//! Time-slice routing and cross-slice likelihoods. Numeric inputs are plain
//! buffers. The R bridge owns input validation, public objects and conditions.
use super::*;

const MODEL_BUFFERS: usize = 22;

struct Buffers<'a> {
    remaining: &'a [&'a [f64]],
}

impl<'a> Buffers<'a> {
    fn take(&mut self, count: usize) -> Result<&'a [&'a [f64]], String> {
        if self.remaining.len() < count {
            return Err("incomplete stratified buffers".into());
        }
        let (head, tail) = self.remaining.split_at(count);
        self.remaining = tail;
        Ok(head)
    }

    fn next(&mut self) -> Result<&'a [f64], String> {
        Ok(self.take(1)?[0])
    }

    fn finish(self) -> Result<(), String> {
        if self.remaining.is_empty() {
            Ok(())
        } else {
            Err("extra stratified buffers".into())
        }
    }
}

fn count(value: f64) -> Result<usize, String> {
    index(value, i32::MAX as usize)
}

fn indices(values: &[f64], limit: usize) -> Result<Vec<usize>, String> {
    values.iter().map(|&v| index(v, limit)).collect()
}

#[derive(Clone, Copy)]
struct Route {
    piece: usize,
    tip: usize,
}

impl Route {
    fn decode(buffer: &[f64]) -> Result<Option<Self>, String> {
        match buffer {
            [] => Ok(None),
            [piece, tip] => Ok(Some(Self {
                piece: count(*piece)?,
                tip: count(*tip)?,
            })),
            _ => Err("invalid stratified route".into()),
        }
    }
}

struct Branch {
    source: usize,
    table_row: usize,
    original_length: f64,
    down_length: Option<f64>,
    up_length: Option<f64>,
    older: Option<Route>,
}

struct Subtree {
    tree: CompiledModelContext,
    sources: Vec<usize>,
    table_rows: Vec<usize>,
    master_nodes: Vec<usize>,
    root_rows: Vec<bool>,
    older: Option<Route>,
    up_fixes: Vec<(usize, Vec<f64>)>,
}

enum Piece {
    Branch(Branch),
    Subtree(Box<Subtree>),
}

struct Layer {
    columns: Vec<usize>,
    geo_mask: Vec<bool>,
    q: Option<CompiledQTemplate>,
    clado: Option<bgb_clado::CompiledCladoTemplate>,
    pieces: Vec<Piece>,
}

struct TipBoundary {
    top: Vec<f64>,
    bottom: Vec<f64>,
}

struct Scratch {
    cond_table: Vec<f64>,
    down_table: Vec<f64>,
    all_cond: Vec<f64>,
    tips: Vec<f64>,
    next_tips: Vec<f64>,
}

pub struct CompiledStratifiedContext {
    n: usize,
    tip_count: usize,
    node_count: usize,
    table_count: usize,
    root: usize,
    null: usize,
    ttl_rows: Vec<usize>,
    down_rows: Vec<usize>,
    tip_routes: Vec<(usize, Route)>,
    layers: Vec<Layer>,
    root_prior: Vec<f64>,
    pool: Arc<rayon::ThreadPool>,
    scratch: Scratch,
}

struct Options {
    ancestors: bool,
    full: bool,
    keep_table: bool,
    ttl_from_table: bool,
    rootedge: bool,
}

impl Options {
    fn decode(a: &[f64]) -> Result<Self, String> {
        if a.len() != 5 {
            return Err("invalid stratified options".into());
        }
        Ok(Self {
            ancestors: a[0] != 0.0,
            full: a[1] != 0.0,
            keep_table: a[2] != 0.0 || a[3] != 0.0,
            ttl_from_table: a[3] != 0.0,
            rootedge: a[4] != 0.0,
        })
    }

    fn append_condlikes(&self) -> bool {
        self.full || !self.ttl_from_table
    }
}

impl CompiledStratifiedContext {
    pub fn new(a: &[&[f64]], threads: usize) -> Result<Self, String> {
        let mut input = Buffers { remaining: a };
        let [n, tips, nodes, rows, root, layers, null] = input.next()? else {
            return Err("invalid stratified context header".into());
        };
        let (n, tips, nodes, rows, layers) = (
            count(*n)?,
            count(*tips)?,
            count(*nodes)?,
            count(*rows)?,
            count(*layers)?,
        );
        let root = index(*root, nodes)?;
        let null = index(*null, 2)?;
        if n == 0 || tips < 2 || nodes != 2 * tips - 1 || rows < tips || layers < 2 {
            return Err("invalid master tree dimensions".into());
        }
        let ttl_rows = indices(input.next()?, rows)?;
        let down_rows = indices(input.next()?, rows)?;
        let encoded_tips = input.next()?;
        if encoded_tips.len() != tips * 3 {
            return Err("invalid master tip routes".into());
        }
        let tip_routes = encoded_tips
            .chunks_exact(3)
            .map(|row| {
                Ok((
                    index(row[0], layers)?,
                    Route {
                        piece: count(row[1])?,
                        tip: count(row[2])?,
                    },
                ))
            })
            .collect::<Result<Vec<_>, String>>()?;
        let pool = Arc::new(
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads.max(1))
                .build()
                .map_err(|e| e.to_string())?,
        );
        let mut specs = Vec::with_capacity(layers);
        let mut current_tips = tips;
        for _ in 0..layers {
            let [piece_count, area_count] = input.next()? else {
                return Err("invalid stratum header".into());
            };
            let pieces = count(*piece_count)?;
            let areas = count(*area_count)?;
            let columns = indices(input.next()?, n)?;
            let geo_mask: Vec<_> = input.next()?.iter().map(|&value| value != 0.0).collect();
            if columns.is_empty() || columns.len() < null || geo_mask.is_empty() || pieces == 0 {
                return Err("empty stratified layer".into());
            }
            let q_flat = input.next()?;
            let q_lens = input.next()?;
            let clado_flat = input.next()?;
            let clado_lens = input.next()?;
            let q = if q_lens.is_empty() {
                None
            } else {
                Some(CompiledQTemplate::new(
                    &unpack_states(q_flat, q_lens)?,
                    areas,
                )?)
            };
            let clado = if clado_lens.is_empty() {
                None
            } else {
                Some(compile_clado(clado_flat, clado_lens, areas)?)
            };
            let mut layer_pieces = Vec::with_capacity(pieces);
            for _ in 0..pieces {
                let meta = input.next()?;
                let piece = match meta {
                    [0.0, slot, row, original, down, down_len, up, up_len] => {
                        Piece::Branch(Branch {
                            source: index(*slot, current_tips)?,
                            table_row: index(*row, rows)?,
                            original_length: *original,
                            down_length: (*down != 0.0).then_some(*down_len),
                            up_length: (*up != 0.0).then_some(*up_len),
                            older: Route::decode(input.next()?)?,
                        })
                    }
                    [1.0] => {
                        let sources = indices(input.next()?, current_tips)?;
                        let table_rows = indices(input.next()?, rows)?;
                        let master_nodes = indices(input.next()?, nodes)?;
                        let root_rows: Vec<_> = input.next()?.iter().map(|&v| v != 0.0).collect();
                        let older = Route::decode(input.next()?)?;
                        let up_nodes = input.next()?;
                        let up_values = input.next()?;
                        let tree = CompiledModelContext::new_with_pool(
                            input.take(13)?,
                            threads,
                            Some(pool.clone()),
                        )?;
                        if tree.n != columns.len()
                            || sources.len() != tree.tips
                            || table_rows.len() != tree.nodes
                            || master_nodes.len() != tree.nodes
                            || root_rows.len() != tree.nodes
                            || up_values.len() != up_nodes.len() * tree.n
                        {
                            return Err("invalid subtree routes or fixed likes".into());
                        }
                        let up_fixes = up_nodes
                            .iter()
                            .enumerate()
                            .map(|(k, &node)| {
                                let values = (0..tree.n)
                                    .map(|j| up_values[k + up_nodes.len() * j])
                                    .collect();
                                Ok((index(node, tree.nodes)?, values))
                            })
                            .collect::<Result<_, String>>()?;
                        Piece::Subtree(Box::new(Subtree {
                            tree,
                            sources,
                            table_rows,
                            master_nodes,
                            root_rows,
                            older,
                            up_fixes,
                        }))
                    }
                    _ => return Err("invalid stratified piece".into()),
                };
                layer_pieces.push(piece);
            }
            specs.push(Layer {
                columns,
                geo_mask,
                q,
                clado,
                pieces: layer_pieces,
            });
            current_tips = pieces;
        }
        input.finish()?;
        let mut root_prior = Vec::new();
        if let Some(oldest) = specs.last() {
            for piece in &oldest.pieces {
                if let Piece::Subtree(sub) = piece {
                    if sub.older.is_none()
                        && sub.master_nodes.get(sub.tree.root) == Some(&root)
                        && !sub.tree.prior.is_empty()
                    {
                        root_prior.resize(n, 0.0);
                        write_columns(&mut root_prior, &sub.tree.prior, &oldest.columns);
                        break;
                    }
                }
            }
        }
        Ok(Self {
            n,
            tip_count: tips,
            node_count: nodes,
            table_count: rows,
            root,
            null,
            ttl_rows,
            down_rows,
            tip_routes,
            layers: specs,
            root_prior,
            pool,
            scratch: Scratch {
                cond_table: vec![0.0; rows * n],
                down_table: Vec::new(),
                all_cond: Vec::new(),
                tips: Vec::new(),
                next_tips: Vec::new(),
            },
        })
    }

    fn validate_contract_context(&self, contract: u32) -> Result<(), String> {
        let strict_contract = contract != 0 || !self.root_prior.is_empty();
        if self.null != 1 || self.tip_routes.iter().any(|&(layer, _)| layer != 0) {
            return Err(
                "UnsupportedModelContractError: stratified probability contracts require include_null_range=true and contemporaneous present-day tips".into(),
            );
        }
        let oldest_idx = self.layers.len() - 1;
        for (i, layer) in self.layers.iter().enumerate() {
            if layer.q.is_none()
                || layer.clado.is_none()
                || layer.geo_mask.len() != self.n
                || !layer.geo_mask[0]
                || layer.columns.len() <= 1
            {
                return Err(
                    "UnsupportedModelContractError: stratified probability contracts exclude trait models and require nonempty geographic strata with an allowed empty state".into(),
                );
            }
            for piece in &layer.pieces {
                match piece {
                    Piece::Branch(b) => {
                        let valid_up = if strict_contract {
                            b.up_length == Some(b.original_length)
                        } else {
                            b.up_length.is_none_or(|u| u == b.original_length)
                        };
                        if b.down_length != Some(b.original_length)
                            || !valid_up
                            || !b.original_length.is_finite()
                            || b.original_length <= 0.0
                        {
                            return Err(
                                "UnsupportedModelContractError: stratified probability contracts reject fossil-hook or truncated branches".into(),
                            );
                        }
                    }
                    Piece::Subtree(s) => {
                        let master_is_root =
                            s.older.is_none() && s.master_nodes.get(s.tree.root) == Some(&self.root);
                        if !s.up_fixes.is_empty()
                            || !s.tree.fixed.iter().all(Vec::is_empty)
                            || s.tree.trait_len != 0
                            || s.tree.kind != 1
                            || s.tree.method != 3
                            || s.tree.steps.iter().any(|st| s.tree.hook(st))
                            || (i == oldest_idx && (!master_is_root || s.tree.root_edge.is_some()))
                        {
                            return Err(
                                "UnsupportedModelContractError: stratified probability contracts reject fixed internal nodes, fossil hooks and root edges".into(),
                            );
                        }
                    }
                }
            }
        }
        for i in 0..oldest_idx {
            let young = &self.layers[i].geo_mask;
            let old = &self.layers[i + 1].geo_mask;
            if old.iter().zip(young).any(|(&allowed_old, &allowed_young)| allowed_old && !allowed_young)
            {
                return Err(
                    "UnsupportedModelContractError: probability contracts require monotone state-space emergence across time strata".into(),
                );
            }
        }
        Ok(())
    }

    pub fn evaluate_contract(
        &mut self,
        a: &[&[f64]],
        mode: PrecisionMode,
        fused: bool,
        contract: u32,
        replacement_tips: &[f64],
    ) -> Result<Vec<Vec<f64>>, String> {
        if contract > 2 {
            return Err("unknown likelihood contract".into());
        }
        self.validate_contract_context(contract)?;
        if a.len() < 2 {
            return Err("incomplete stratified buffers".into());
        }
        let candidate_tips = if replacement_tips.is_empty() {
            a[1]
        } else {
            replacement_tips
        };
        if candidate_tips.len() != self.tip_count * self.n
            || candidate_tips.iter().any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err("invalid tip likelihood vector".into());
        }
        let present_mask = &self.layers[0].geo_mask;
        if matches!(contract, 1 | 2) {
            for i in 0..self.tip_count {
                if (0..self.n).any(|s| candidate_tips[i + self.tip_count * s] > 0.0 && !present_mask[s])
                {
                    return Err(
                        "tip likelihoods assign positive mass to a state disallowed in the present stratum".into(),
                    );
                }
                if !(0..self.n).any(|s| candidate_tips[i + self.tip_count * s] > 0.0)
                    || (contract == 1 && candidate_tips[i] != 0.0)
                {
                    return Err("tip likelihoods violate probability contract".into());
                }
            }
        }
        let oldest_mask = &self.layers.last().unwrap().geo_mask;
        if self.root_prior.is_empty() {
            if contract != 0 {
                return Err("missing root likelihood weights for probability contract".into());
            }
        } else if self.root_prior.len() != self.n
            || self.root_prior.iter().any(|x| !x.is_finite() || *x < 0.0)
            || self.root_prior.iter().sum::<f64>() <= 0.0
            || self
                .root_prior
                .iter()
                .zip(oldest_mask)
                .any(|(&w, &ok)| w > 0.0 && !ok)
            || (contract == 1 && self.root_prior[0] != 0.0)
        {
            return Err("invalid root likelihood weights".into());
        }
        self.evaluate_inner(
            a,
            Some(candidate_tips),
            mode,
            fused,
            if contract == 0 { 3 } else { contract },
        )
    }

    pub fn evaluate(
        &mut self,
        a: &[&[f64]],
        mode: PrecisionMode,
        fused: bool,
    ) -> Result<Vec<Vec<f64>>, String> {
        self.evaluate_inner(a, None, mode, fused, 0)
    }

    fn run_downpass_only(
        &mut self,
        layer_inputs: &[&[&[f64]]],
        tip_input: &[f64],
        mode: PrecisionMode,
        fused: bool,
        internal_contract: u32,
    ) -> Result<f64, String> {
        let options = Options {
            ancestors: false,
            full: false,
            keep_table: true,
            ttl_from_table: true,
            rootedge: false,
        };
        let n = self.n;
        let scratch = &mut self.scratch;
        scratch.cond_table.fill(0.0);
        scratch.all_cond.clear();
        scratch.tips.resize(tip_input.len(), 0.0);
        for tip in 0..self.tip_count {
            let row = row_mut(&mut scratch.tips, tip, n);
            for j in 0..n {
                row[j] = tip_input[tip + self.tip_count * j];
            }
            row_mut(&mut scratch.cond_table, tip, n).copy_from_slice(row);
            normalize(row);
        }
        for (layer, &down) in self.layers.iter_mut().zip(layer_inputs) {
            let width = layer.columns.len();
            let propagation = build_propagator(width, layer.q.as_ref(), down, mode, fused)?;
            let mut events = build_clado(width, self.null, 1, 3, layer.clado.as_ref(), down, mode)?;
            events.probability_contract = internal_contract != 0;
            events.complete_empty = internal_contract == 2;
            let layer_down = self.pool.install(|| {
                layer
                    .pieces
                    .par_iter_mut()
                    .map(|piece| {
                        piece.downpass(
                            n,
                            &layer.columns,
                            &layer.geo_mask,
                            &scratch.tips,
                            &propagation,
                            &events,
                            down[0][1] != 0.0,
                            down[0][3] as usize,
                            &options,
                        )
                    })
                    .collect::<Result<Vec<_>, String>>()
            })?;
            scratch.next_tips.resize(layer.pieces.len() * n, 0.0);
            for (j, (piece, result)) in layer.pieces.iter().zip(layer_down).enumerate() {
                piece.commit_down(&layer.columns, result, j, n, scratch, &options)?;
            }
            std::mem::swap(&mut scratch.tips, &mut scratch.next_tips);
        }
        if self.ttl_rows.len() != self.node_count {
            return Err("master likelihood row count".into());
        }
        log_row_sums(
            self.ttl_rows
                .iter()
                .map(|&i| row(&scratch.cond_table, i, n)),
        )
    }

    fn evaluate_inner(
        &mut self,
        a: &[&[f64]],
        tip_override: Option<&[f64]>,
        mode: PrecisionMode,
        fused: bool,
        internal_contract: u32,
    ) -> Result<Vec<Vec<f64>>, String> {
        let mut input = Buffers { remaining: a };
        let options = Options::decode(input.next()?)?;
        let raw_tip_input = input.next()?;
        let tip_input = tip_override.unwrap_or(raw_tip_input);
        if tip_input.len() != self.tip_count * self.n
            || options.ancestors && !options.ttl_from_table
            || (internal_contract != 0 && (!options.keep_table || !options.ttl_from_table))
        {
            return Err("unsupported stratified tip dimensions or ancestor options".into());
        }

        let mut layer_pairs = Vec::with_capacity(self.layers.len());
        for _ in 0..self.layers.len() {
            let [separate] = input.next()? else {
                return Err("missing uppass switch".into());
            };
            let down = input.take(MODEL_BUFFERS)?;
            let up = if *separate != 0.0 {
                input.take(MODEL_BUFFERS)?
            } else {
                down
            };
            if down[0].len() != 4 || up[0].len() != 4 {
                return Err("invalid stratum model flags".into());
            }
            layer_pairs.push((*separate != 0.0, down, up));
        }
        input.finish()?;

        let log_survival = if internal_contract == 1 {
            let present_mask = &self.layers[0].geo_mask;
            let mut survival_tips = vec![0.0; self.tip_count * self.n];
            for s in self.null..self.n {
                if present_mask[s] {
                    for tip in 0..self.tip_count {
                        survival_tips[tip + self.tip_count * s] = 1.0;
                    }
                }
            }
            let down_only: Vec<_> = layer_pairs.iter().map(|&(_, down, _)| down).collect();
            let value =
                self.run_downpass_only(&down_only, &survival_tips, mode, fused, internal_contract)?;
            if !value.is_finite() {
                return Err("zero or invalid tip-survival probability".into());
            }
            value
        } else {
            f64::NAN
        };

        let n = self.n;
        let scratch = &mut self.scratch;
        scratch.cond_table.fill(0.0);
        scratch.all_cond.clear();
        if options.ancestors {
            scratch.down_table.resize(self.table_count * n, 0.0);
            scratch.down_table.fill(0.0);
        }
        scratch.tips.resize(tip_input.len(), 0.0);
        for tip in 0..self.tip_count {
            let row = row_mut(&mut scratch.tips, tip, n);
            for j in 0..n {
                row[j] = tip_input[tip + self.tip_count * j];
            }
            if options.keep_table {
                row_mut(&mut scratch.cond_table, tip, n).copy_from_slice(row);
            }
            normalize(row);
            if options.append_condlikes() {
                scratch.all_cond.extend_from_slice(row);
            }
        }
        let mut upper_models = Vec::with_capacity(self.layers.len());
        for (layer, (separate, down, up)) in self.layers.iter_mut().zip(layer_pairs) {
            let width = layer.columns.len();
            let propagation = build_propagator(width, layer.q.as_ref(), down, mode, fused)?;
            let mut events = build_clado(width, self.null, 1, 3, layer.clado.as_ref(), down, mode)?;
            events.probability_contract = internal_contract != 0;
            events.complete_empty = internal_contract == 2;
            let layer_down = self.pool.install(|| {
                layer
                    .pieces
                    .par_iter_mut()
                    .map(|piece| {
                        piece.downpass(
                            n,
                            &layer.columns,
                            &layer.geo_mask,
                            &scratch.tips,
                            &propagation,
                            &events,
                            down[0][1] != 0.0,
                            down[0][3] as usize,
                            &options,
                        )
                    })
                    .collect::<Result<Vec<_>, String>>()
            })?;
            scratch.next_tips.resize(layer.pieces.len() * n, 0.0);
            for (j, (piece, result)) in layer.pieces.iter().zip(layer_down).enumerate() {
                piece.commit_down(&layer.columns, result, j, n, scratch, &options)?;
            }
            std::mem::swap(&mut scratch.tips, &mut scratch.next_tips);
            if options.ancestors {
                let same_q = down[0][..2] == up[0][..2]
                    && [1, 2, 3, 10, 11, 12, 13].iter().all(|&j| down[j] == up[j]);
                let up_prop = if same_q {
                    propagation
                } else {
                    build_propagator(width, layer.q.as_ref(), up, mode, fused)?
                };
                let mut up_events = if !separate {
                    events
                } else {
                    let mut built =
                        build_clado(width, self.null, 1, 3, layer.clado.as_ref(), up, mode)?;
                    built.probability_contract = internal_contract != 0;
                    built.complete_empty = internal_contract == 2;
                    built
                };
                // The historical stratified hook constructs diagonal events in
                // geographic state space even for trait-expanded probabilities.
                up_events.hook_states = layer.geo_mask.iter().filter(|&&x| x).count() - self.null;
                upper_models.push((up_prop, up_events));
            }
        }

        let raw_likelihood = if options.ttl_from_table {
            if self.ttl_rows.len() != self.node_count {
                return Err("master likelihood row count".into());
            }
            if internal_contract != 0
                && self
                    .ttl_rows
                    .iter()
                    .any(|&i| sum(row(&scratch.cond_table, i, n)) < 0.0)
            {
                return Err("negative likelihood requires original R warning semantics".into());
            }
            log_row_sums(
                self.ttl_rows
                    .iter()
                    .map(|&i| row(&scratch.cond_table, i, n)),
            )?
        } else {
            let count = scratch
                .all_cond
                .chunks_exact(n)
                .filter(|r| r.iter().any(|&x| x != 0.0))
                .count();
            log_row_sums(
                scratch
                    .all_cond
                    .chunks_exact(n)
                    .filter(|r| r.iter().any(|&x| x != 0.0))
                    .take(count.saturating_sub(usize::from(!options.rootedge))),
            )?
        };
        let likelihood = raw_likelihood - if internal_contract == 1 { log_survival } else { 0.0 };
        if internal_contract != 0 && options.ancestors && !likelihood.is_finite() {
            return Err("ancestral probabilities are undefined for zero-probability data".into());
        }
        let emit_contract_meta = matches!(internal_contract, 1 | 2);
        if !options.full && !options.ancestors {
            return Ok(vec![if emit_contract_meta {
                vec![likelihood, log_survival]
            } else {
                vec![likelihood]
            }]);
        }

        let all_condlikes: Vec<_> = scratch
            .all_cond
            .chunks_exact(n)
            .filter(|r| r.iter().any(|&v| v != 0.0))
            .flatten()
            .copied()
            .collect();
        let mut out = vec![
            vec![likelihood, (all_condlikes.len() / n) as f64],
            column_major(&all_condlikes, all_condlikes.len() / n, n),
            column_major(&scratch.cond_table, self.table_count, n),
        ];
        if options.ancestors {
            if self.down_rows.len() + 1 != self.node_count {
                return Err("master downpass row count".into());
            }
            let mut top_down = Vec::with_capacity(self.node_count * n);
            for &i in &self.ttl_rows {
                top_down.extend(normalized(row(&scratch.cond_table, i, n)));
            }
            let mut bottom_down = Vec::with_capacity(self.node_count * n);
            for (i, &r) in self.down_rows.iter().enumerate() {
                if i == self.tip_count {
                    bottom_down.extend(std::iter::repeat_n(missing(), n));
                }
                bottom_down.extend(normalized(row(&scratch.down_table, r, n)));
            }
            if bottom_down.len() != self.node_count * n {
                return Err("master downpass root row missing".into());
            }
            let down_table = column_major(&scratch.down_table, self.table_count, n);
            let (top_up, bottom_up) =
                self.uppass(upper_models, &bottom_down, internal_contract != 0)?;
            let top_marginal = marginal(&top_down, &top_up, self.root, n)?;
            let bottom_marginal = marginal(&bottom_down, &bottom_up, self.root, n)?;
            out.push(down_table);
            for values in [
                &top_up,
                &bottom_up,
                &bottom_marginal,
                &top_marginal,
                &bottom_down,
                &top_down,
            ] {
                out.push(column_major(values, self.node_count, n));
            }
        }
        if emit_contract_meta {
            out[0].push(log_survival);
        }
        Ok(out)
    }

    fn uppass(
        &mut self,
        models: Vec<(Propagator, Clado)>,
        down_bottom: &[f64],
        probability_contract: bool,
    ) -> Result<(Vec<f64>, Vec<f64>), String> {
        let n = self.n;
        let mut top = vec![0.0; self.node_count * n];
        let mut bottom = vec![0.0; self.node_count * n];
        if probability_contract && !self.root_prior.is_empty() {
            row_mut(&mut top, self.root, n).copy_from_slice(&self.root_prior);
        } else {
            row_mut(&mut top, self.root, n).fill(1.0 / n as f64);
        }
        let mut routes: Vec<Vec<TipBoundary>> =
            (0..self.layers.len()).map(|_| Vec::new()).collect();
        let mut last_bottom = None;
        let oldest = self.layers.len() - 1;
        for (i, (layer, (p, c))) in self.layers.iter_mut().zip(models).enumerate().rev() {
            let (current, older) = routes.split_at_mut(i + 1);
            for piece in &mut layer.pieces {
                let previous = match piece {
                    Piece::Branch(branch) => branch.older,
                    Piece::Subtree(tree) => tree.older,
                };
                let prev = previous
                    .map(|route| boundary(older.first().map_or(&[][..], Vec::as_slice), route, n))
                    .transpose()?;
                let piece_boundary = match piece {
                    Piece::Branch(branch) => {
                        let (start, bot) =
                            prev.ok_or("stratified branch missing older boundary")?;
                        last_bottom = Some(bot.to_vec());
                        let output = if let Some(length) = branch.up_length {
                            propagate_boundary(
                                &p,
                                length,
                                start,
                                &layer.columns,
                                &layer.geo_mask,
                                self.null,
                                c.probability_contract,
                            )?
                        } else {
                            start.to_vec()
                        };
                        TipBoundary {
                            top: output,
                            bottom: bot.to_vec(),
                        }
                    }
                    Piece::Subtree(piece) => {
                        let tree = &mut piece.tree;
                        let root_master = piece.master_nodes[tree.root];
                        let (seed_top, seed_bottom) = if i == oldest {
                            (
                                gather(row(&top, root_master, n), &layer.columns),
                                vec![0.0; tree.n],
                            )
                        } else if let Some(length) = tree.root_edge {
                            let (start, bot) =
                                prev.ok_or("stratified subtree missing older boundary")?;
                            last_bottom = Some(bot.to_vec());
                            let output = propagate_boundary(
                                &p,
                                length,
                                start,
                                &layer.columns,
                                &layer.geo_mask,
                                self.null,
                                c.probability_contract,
                            )?;
                            (gather(&output, &layer.columns), gather(bot, &layer.columns))
                        } else {
                            // Upstream reads the existing binding for this branch.
                            let bot = last_bottom.as_ref().ok_or("unbound prior branch bottom")?;
                            (
                                gather(row(&top, root_master, n), &layer.columns),
                                gather(bot, &layer.columns),
                            )
                        };
                        tree.bottom.resize(tree.nodes * tree.n, 0.0);
                        for (j, &master) in piece.master_nodes.iter().enumerate() {
                            for (k, &col) in layer.columns.iter().enumerate() {
                                tree.bottom[j * tree.n + k] = down_bottom[master * n + col];
                            }
                        }
                        // Uppass fixes differ from downpass fixes for trait models.
                        let mut fixes = vec![Vec::new(); tree.nodes];
                        for (node, values) in &piece.up_fixes {
                            fixes[*node] = values.clone();
                        }
                        std::mem::swap(&mut tree.fixed, &mut fixes);
                        let evaluated = tree.uppass(&p, &c, Some((&seed_top, &seed_bottom)));
                        std::mem::swap(&mut tree.fixed, &mut fixes);
                        evaluated?;
                        for (j, &master) in piece.master_nodes.iter().enumerate() {
                            write_columns(
                                row_mut(&mut top, master, n),
                                tree.row(&tree.up, j),
                                &layer.columns,
                            );
                            write_columns(
                                row_mut(&mut bottom, master, n),
                                tree.row(&tree.upbot, j),
                                &layer.columns,
                            );
                        }
                        let mut out = TipBoundary {
                            top: vec![0.0; tree.tips * n],
                            bottom: vec![0.0; tree.tips * n],
                        };
                        for tip in 0..tree.tips {
                            write_columns(
                                row_mut(&mut out.top, tip, n),
                                tree.row(&tree.up, tip),
                                &layer.columns,
                            );
                            write_columns(
                                row_mut(&mut out.bottom, tip, n),
                                tree.row(&tree.upbot, tip),
                                &layer.columns,
                            );
                        }
                        out
                    }
                };
                if piece_boundary.top.iter().any(|x| x.is_nan()) {
                    return Err("stratified boundary NaN".into());
                }
                current[i].push(piece_boundary);
            }
        }
        for (j, &(layer, route)) in self.tip_routes.iter().enumerate() {
            let (start, _) = boundary(&routes[layer], route, n)?;
            row_mut(&mut top, j, n).copy_from_slice(start);
        }
        Ok((top, bottom))
    }
}

impl Piece {
    fn downpass(
        &mut self,
        n: usize,
        columns: &[usize],
        geo_mask: &[bool],
        tips: &[f64],
        p: &Propagator,
        c: &Clado,
        sparse: bool,
        trait_len: usize,
        options: &Options,
    ) -> Result<(Vec<f64>, Vec<f64>), String> {
        match self {
            Self::Branch(branch) => {
                let start = row(tips, branch.source, n);
                let time = if options.keep_table {
                    branch.down_length
                } else {
                    Some(branch.original_length)
                };
                let result = if let Some(t) = time {
                    let local = p.apply(t, &gather(start, columns), false)?;
                    let mut result = vec![0.0; n];
                    write_columns(&mut result, &local, columns);
                    mask_geography(&mut result, geo_mask);
                    result
                } else {
                    start.to_vec()
                };
                Ok((result, Vec::new()))
            }
            Self::Subtree(piece) => {
                let tree = &mut piece.tree;
                tree.trait_len = trait_len;
                for (i, &src) in piece.sources.iter().enumerate() {
                    for (j, &col) in columns.iter().enumerate() {
                        tree.tip_likes[i + tree.tips * j] = tips[src * n + col];
                    }
                }
                tree.downpass(p, c, sparse, options.ancestors)
            }
        }
    }

    fn commit_down(
        &self,
        columns: &[usize],
        result: (Vec<f64>, Vec<f64>),
        next_slot: usize,
        n: usize,
        scratch: &mut Scratch,
        options: &Options,
    ) -> Result<(), String> {
        let (conditional, root_bottom) = result;
        let next = row_mut(&mut scratch.next_tips, next_slot, n);
        next.fill(0.0);
        match self {
            Self::Branch(branch) => {
                next.copy_from_slice(&conditional);
                if options.keep_table {
                    row_mut(&mut scratch.cond_table, branch.table_row, n).copy_from_slice(next);
                }
                if options.ancestors {
                    row_mut(&mut scratch.down_table, branch.table_row, n)
                        .copy_from_slice(&normalized(next));
                }
            }
            Self::Subtree(piece) => {
                let tree = &piece.tree;
                let root_cond = if conditional.is_empty() {
                    tree.row(&tree.cond, tree.root)
                } else {
                    &conditional
                };
                write_columns(next, root_cond, columns);
                for (j, &table_row) in piece.table_rows.iter().enumerate() {
                    let row_cond = tree.row(&tree.cond, j);
                    if options.keep_table {
                        write_columns(
                            row_mut(&mut scratch.cond_table, table_row, n),
                            row_cond,
                            columns,
                        );
                    }
                    if options.ancestors {
                        let down = if piece.root_rows[j] {
                            &root_bottom
                        } else {
                            tree.row(&tree.bottom, j)
                        };
                        write_columns(
                            row_mut(&mut scratch.down_table, table_row, n),
                            down,
                            columns,
                        );
                    }
                    if j >= tree.tips && options.append_condlikes() {
                        append_scattered(&mut scratch.all_cond, row_cond, columns, n);
                    }
                }
                if !conditional.is_empty() && options.append_condlikes() {
                    append_scattered(&mut scratch.all_cond, &conditional, columns, n);
                }
            }
        }
        if next.iter().any(|x| x.is_nan()) {
            return Err("stratified downpass NaN".into());
        }
        Ok(())
    }
}

fn row(values: &[f64], i: usize, width: usize) -> &[f64] {
    &values[i * width..(i + 1) * width]
}

fn row_mut(values: &mut [f64], i: usize, width: usize) -> &mut [f64] {
    &mut values[i * width..(i + 1) * width]
}

fn gather(values: &[f64], columns: &[usize]) -> Vec<f64> {
    columns.iter().map(|&col| values[col]).collect()
}

fn write_columns(target: &mut [f64], values: &[f64], columns: &[usize]) {
    for (&col, &value) in columns.iter().zip(values) {
        target[col] = value;
    }
}

fn append_scattered(target: &mut Vec<f64>, values: &[f64], columns: &[usize], n: usize) {
    let start = target.len();
    target.resize(start + n, 0.0);
    write_columns(&mut target[start..], values, columns);
}

fn normalize(values: &mut [f64]) {
    let total = sum(values);
    if total <= 0.0 || !total.is_finite() {
        values.fill(0.0);
        return;
    }
    for value in values {
        *value /= total;
    }
}

fn normalized(values: &[f64]) -> Vec<f64> {
    let mut result = values.to_vec();
    normalize(&mut result);
    result
}

fn mask_geography(values: &mut [f64], mask: &[bool]) {
    for (value, &allowed) in values.iter_mut().zip(mask.iter().cycle()) {
        if !allowed {
            *value = 0.0;
        }
    }
}

fn boundary(pieces: &[TipBoundary], route: Route, n: usize) -> Result<(&[f64], &[f64]), String> {
    let piece = pieces.get(route.piece).ok_or("invalid older piece")?;
    if route.tip >= piece.top.len() / n {
        return Err("invalid boundary tip".into());
    }
    Ok((
        row(&piece.top, route.tip, n),
        row(&piece.bottom, route.tip, n),
    ))
}

fn propagate_boundary(
    p: &Propagator,
    length: f64,
    start: &[f64],
    columns: &[usize],
    geo_mask: &[bool],
    null: usize,
    probability_contract: bool,
) -> Result<Vec<f64>, String> {
    let mut local = p.apply(length, &gather(start, columns), true)?;
    if !probability_contract && null != 0 {
        local[0] = 0.0;
    }
    normalize(&mut local);
    let mut result = vec![0.0; start.len()];
    write_columns(&mut result, &local, columns);
    mask_geography(&mut result, geo_mask);
    normalize(&mut result);
    Ok(result)
}

fn log_row_sums<'a>(rows: impl Iterator<Item = &'a [f64]>) -> Result<f64, String> {
    let logs: Vec<_> = rows.map(|r| sum(r).ln()).collect();
    if logs.iter().any(|x| x.is_nan()) {
        return Err("stratified likelihood NaN requires reference diagnostics".into());
    }
    Ok(sum(&logs))
}

fn marginal(down: &[f64], up: &[f64], root: usize, n: usize) -> Result<Vec<f64>, String> {
    let mut result = Vec::with_capacity(down.len());
    for (node, (d, u)) in down.chunks_exact(n).zip(up.chunks_exact(n)).enumerate() {
        if node == root {
            result.extend(normalized(d));
        } else {
            let start = result.len();
            result.extend(d.iter().zip(u).map(|(a, b)| a * b));
            normalize(&mut result[start..]);
            if result[start..].iter().any(|x| x.is_nan()) {
                return Err("stratified marginal NaN requires reference diagnostics".into());
            }
        }
    }
    Ok(result)
}

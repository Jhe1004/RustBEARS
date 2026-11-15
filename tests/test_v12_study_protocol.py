"""Scientific reporting must retain failed datasets and cluster by dataset."""
import json

import pytest

from scripts.run_v12_v3_study import summarize,verify_prepared_study,verify_completed_dataset


def score(mcc, multiple):
    return dict(support_bin=">0.95",truth_state=0,
        mcc=dict(brier=mcc,log_score=.2,zero_probability=False,probabilities=[.8,.2]),
        multitree=dict(brier=multiple,log_score=.1,zero_probability=False,probabilities=[.9,.1]))


def test_v3_clusters_whole_datasets_and_retains_missing(tmp_path):
    # Ten repeated nodes in dataset 1 must not outweigh one node in dataset 2.
    for index,rows in ((1,[score(.6,.2)]*10),(2,[score(.1,.3)])):
        folder=tmp_path/f'dataset_{index:03d}';folder.mkdir()
        (folder/'v3_scores.json').write_text(json.dumps(dict(paired_clades=rows)))
    records=[dict(dataset=1,status='complete'),dict(dataset=2,status='complete'),
             dict(dataset=3,status='failed_retained',stage='diagnostics_and_import')]
    result=summarize(records,tmp_path,20261115)
    assert result['n_planned']==3 and result['n_missing']==1
    cell=result['support_bins']['all']
    assert cell['n_datasets_with_scored_clades']==2
    assert cell['n_clades']==11
    assert cell['paired_dataset_mean_brier_difference']==pytest.approx(-.1)
    assert cell['dataset_mean_brier']['mcc']==pytest.approx(.35)
    assert cell['brier_missing_result_bounds']['mcc']==pytest.approx([.7/3,2.7/3])
    assert cell['paired_brier_missing_result_bounds']==pytest.approx([-2.2/3,1.8/3])
    assert result['support_bins']['<0.7']['paired_dataset_mean_brier_difference'] is None


def test_v3_zero_probability_is_not_silently_removed(tmp_path):
    folder=tmp_path/'dataset_001';folder.mkdir()
    row=score(.4,2.);row['multitree'].update(log_score=None,zero_probability=True,probabilities=[0.,1.])
    (folder/'v3_scores.json').write_text(json.dumps(dict(paired_clades=[row])))
    result=summarize([dict(dataset=1,status='complete')],tmp_path,20261115)
    cell=result['support_bins']['all']
    assert cell['paired_dataset_bootstrap_ci95'] is None
    assert cell['dataset_scores'][0]['multitree']['infinite_log_score'] is True
    assert cell['dataset_scores'][0]['multitree']['mean_negative_log_probability'] is None
    assert cell['negative_log_probability']['multitree']['infinite'] is True
    assert cell['calibration_bins']['multitree'][0]['sum_observed']==1


def test_v3_cannot_resume_pilot_as_formal(tmp_path):
    job=dict(replicates=1,tips=50,areas=4,sites=2000)
    (tmp_path/'study_plan.json').write_text(json.dumps(dict(mode='pilot',requested_dimensions=job)))
    with pytest.raises(ValueError,match='pilot promotion'):
        verify_prepared_study(tmp_path,job,'locked-budget-sha')


def test_v3_modified_completed_score_is_rejected(tmp_path):
    (tmp_path/'simulation.json').write_text(json.dumps(dict(files={})))
    (tmp_path/'v3_scores.json').write_text('{"tampered":true}')
    with pytest.raises(ValueError,match='score was modified'):
        verify_completed_dataset(tmp_path,dict(score_sha256='original-score-sha'))

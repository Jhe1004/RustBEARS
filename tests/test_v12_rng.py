"""Stable task seeds, declared-batch isolation, and explicit collision handling."""
from dataclasses import replace
import unittest
from unittest.mock import patch

from rustbears.uncertainty.rng import (
    SeedCollisionError, TaskIdentity, canonical_task_payload, derive_seed,
    derive_task_seeds, task_fingerprint,
)


class TaskSeedTests(unittest.TestCase):
    def setUp(self):
        self.task = TaskIdentity("bootstrap_null", "DEC", "1", outer_rep=1, role="simulate")

    def test_cross_language_golden_vector(self):
        self.assertEqual(derive_seed(12345, self.task), 17990196571989855135)
        self.assertEqual(task_fingerprint(12345, self.task),
                         "f9aa047750e7cb9f46518889bfba1ef2c54f399782a2049ea1e4f076404d8806")
        self.assertEqual(canonical_task_payload(12345, self.task),
                         '{"master_seed":12345,"task":{"inner_rep":0,"layer_id":0,"model_id":"DEC","outer_rep":1,"purpose":"bootstrap_null","role":"simulate","step_id":0,"tree_id":"1"}}')

    def test_every_task_dimension_is_separated_and_order_independent(self):
        tasks = [self.task, *(replace(self.task, **{key: value}) for key, value in
                  dict(purpose="bootstrap_alt", model_id="DEC+J", tree_id="2", outer_rep=2,
                       inner_rep=1, role="fit", layer_id=1, step_id=1).items())]
        seeds = derive_task_seeds(12345, tasks)
        self.assertEqual(len(set(seeds.values())), len(tasks))
        self.assertEqual(seeds, derive_task_seeds(12345, reversed(tasks)))
        self.assertNotEqual(seeds[self.task], derive_seed(12346, self.task))

    def test_duplicate_task_and_injected_truncated_collision_are_errors(self):
        with self.assertRaisesRegex(ValueError, "duplicate"):
            derive_task_seeds(12345, [self.task, self.task])
        with patch("rustbears.uncertainty.rng._task_digest", return_value=bytes(32)):
            with self.assertRaises(SeedCollisionError):
                derive_task_seeds(12345, [self.task, replace(self.task, outer_rep=2)])

    def test_integer_limits_and_identifier_validation(self):
        self.assertLessEqual(derive_seed((1 << 64) - 1, self.task), (1 << 64) - 1)
        for value in (True, -1, 1.0, 1 << 64):
            with self.subTest(value=value), self.assertRaises(ValueError):
                derive_seed(value, self.task)
            with self.subTest(rep=value), self.assertRaises(ValueError):
                replace(self.task, outer_rep=value)
        with self.assertRaises(ValueError):
            replace(self.task, model_id="")


if __name__ == "__main__":
    unittest.main()

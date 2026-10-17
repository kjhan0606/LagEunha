"""Run all unit gates with per-test isolation of mutable model configuration.

Build p2_geometry_omp.so before running. Execute on a compute host, not a
login node. --reverse checks the opposite test order; no tests are filtered.
"""
import argparse
import json
import os
from pathlib import Path
import unittest

import moment_ale_evolve as model


class IsolatedResult(unittest.TextTestResult):
    def startTest(self, test):
        self.saved_environment = os.environ.copy()
        self.saved_domain = model.domain
        self.saved_boundary = model.boundary_velocity
        super().startTest(test)

    def stopTest(self, test):
        model.domain = self.saved_domain
        model.boundary_velocity = self.saved_boundary
        os.environ.clear()
        os.environ.update(self.saved_environment)
        super().stopTest(test)


def flatten(suite):
    for item in suite:
        if isinstance(item, unittest.TestSuite):
            yield from flatten(item)
        else:
            yield item


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reverse', action='store_true')
    parser.add_argument('--report', default='verification_suite_report.json')
    args = parser.parse_args()
    suite = unittest.defaultTestLoader.discover(str(Path(__file__).parent), pattern='test_*.py')
    tests = list(flatten(suite))
    if args.reverse:
        tests.reverse()
    result = unittest.TextTestRunner(verbosity=2, resultclass=IsolatedResult).run(unittest.TestSuite(tests))
    report = dict(tests_run=result.testsRun, passed=result.wasSuccessful(),
                  reverse=args.reverse, failures=[str(t) for t, _ in result.failures],
                  errors=[str(t) for t, _ in result.errors],
                  skipped=[(str(t), reason) for t, reason in result.skipped])
    Path(args.report).write_text(json.dumps(report, indent=2) + '\n')
    return 0 if result.wasSuccessful() and not result.skipped else 1


if __name__ == '__main__':
    raise SystemExit(main())

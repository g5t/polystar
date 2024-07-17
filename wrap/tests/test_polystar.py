from unittest import TestCase


class PolystarTestCase(TestCase):

    def test_polystar_importable(self):
        from importlib.util import find_spec
        spec = find_spec('polystar')
        self.assertIsNotNone(spec, "The polystar module could not be found")

        try:
            from importlib import import_module
            import_module('polystar')
        except ModuleNotFoundError as exception:
            self.fail(f"Importing the polystar module failed because it was not found with message\n\t{exception}")
        except ImportError as exception:
            self.fail(f"Importing the polystar module failed with message\n{exception}")



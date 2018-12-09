test-compile-examples
=====================

This test ensures that all models in examples can be compiled.

This currently works by finding any makefile located beneath examples. As such this test will also attempt to build documentation provided with some of the models and hence may fail if a suitable document compiler is not available.
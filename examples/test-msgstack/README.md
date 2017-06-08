test-msgstack
=============

Test if the message stack and TRACE macro work correctly.

This test adds an item to the global message stack, then uses the TRACE macro in
a block. After exiting the block, the message stack should contain only the
first message.

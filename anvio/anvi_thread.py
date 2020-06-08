# coding: utf-8
"""A lightly modified Thread sublcass.  Tracks the return values of Thread's target."""

import anvio

from threading import Thread

__author__ = "Developers of anvi'o (see AUTHORS.txt)"
__copyright__ = "Copyleft 2015-2020, the Meren Lab (http://merenlab.org/)"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = anvio.__version__
__maintainer__ = "Ryan Moore"
__email__ = "moorer@udel.edu"


class AnviThread(Thread):
    """Like regular Thread, except that it tracks the return value of the `target` function call.

    You can then use this return value so the caller can decide what to do.  For example, if you're running a bash
    command, you may want to retry it if the return value is > 0, (i.e., indicating failure).
    """

    def __init__(self, group=None, target=None, name=None, args=(), kwargs=None, *, daemon=None):

        super().__init__(group, target, name, args, kwargs, daemon=daemon)

        # Instance variable to track return value of target function.
        self.target_return_value = None

    def run(self):
        """Method representing the thread's activity.

        You may override this method in a subclass. The standard run() method
        invokes the callable object passed to the object's constructor as the
        target argument, if any, with sequential and keyword arguments taken
        from the args and kwargs arguments, respectively.

        """
        try:
            if self._target:
                # Save the output of the target function in an instance variable for later use.
                self.target_return_value = self._target(*self._args, **self._kwargs)
        finally:
            # Avoid a refcycle if the thread is running a function with
            # an argument that has a member that points to the thread.
            del self._target, self._args, self._kwargs

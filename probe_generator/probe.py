"""Shared objects for probe classes.

"""
class InvalidStatement(Exception):
    """Raised when a probe statement cannot be parsed.

    """


class NonFatalError(Exception):
    """Superclass of exceptions which signal that, although the sequence of the
    probe cannot be determined, processing of the rest of the probes should
    continue.

    """

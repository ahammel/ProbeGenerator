"""Shared exceptions for probe_generator.

"""
class NonFatalError(Exception):
    """Superclass of exceptions which signal that, although the sequence of the
    probe cannot be determined, processing of the rest of the probes should
    continue.

    """

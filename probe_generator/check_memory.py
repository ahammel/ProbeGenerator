"""Find the total memory on a linux system.

"""
def total_ram():
    """Return the total install RAM on the system in kB.

    Raises an Error if the file '/proc/meminfo' is not available or if it
    doesn't contain a 'MemTotal' entry (probably because it's being run on a
    non-Linux machine).

    """
    try:
        with open('/proc/meminfo') as handle:
            for line in handle:
                if line.startswith("MemTotal:"):
                    total_memory = int(line.split()[1])
                    break
            else:
                raise Error("'MemTotal:' not in /proc/meminfo")
    except IOError as error:
        raise Error(str(error))
    return total_memory

class Error(Exception):
    """Generic error class for check_memory module.

    """

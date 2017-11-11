def deduplicate_tuple(t):
    """remove duplicate entries from the given tuple (starting from the left)"""
    seen_entries = []
    result = ()
    for entry in t:
        if entry in seen_entries:
            continue
        seen_entries.append(entry)
        result += (entry,)
    return result


def has_duplicates(t):
    """Return True, if given tuple t has duplicate entries, otherwise return False."""
    for d in _duplicates(t):
        return True
    return False


def _duplicates(t):
    """generator yielding duplicate entries of given tuple (starting from the left)"""
    seen_entries = []
    for entry in t:
        if entry in seen_entries:
            yield entry
        seen_entries.append(entry)

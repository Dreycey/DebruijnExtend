#! usr/bin/python3
"""
This modules defines exceptions caught by the DebruijnExtend tool
"""

# error for unacceptable hash table 
class DebextendError(Exception):
    """Base class Error for DebruijnExtend"""
    pass

# error for unacceptable hash table 
class HashTableError(DebextendError):
    """Base class for other exceptions"""
    pass

class HashTableClusterMismatchError(HashTableError):
    """Base class for other exceptions"""
    pass
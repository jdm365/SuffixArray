# SuffixArray 
Truncated SuffixArray based search index algorithm written in c++ and exposed through cython.

To install
```
make install
```

```python
from suffix_array import SuffixArray

docs = [
    "The quick brown fox jumps over the lazy dog",
    "I am going to the store to buy some milk",
    "Uhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh"
]

## Creates and builds the index
suffix_array = SuffixArray(documents=docs, max_suffix_length=32)


## Query the index. Returns idxs of documents that contain the query substring.
idxs = suffix_array.query("the quick brown fox")
```

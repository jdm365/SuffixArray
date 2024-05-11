# SuffixArray 
Truncated SuffixArray based substring search index algorithm written in c++ and exposed through cython.
Typical query times range from microseconds to milliseconds even on tens or hundreds of millions of items.
**More comprehensive benchmarks to come as package get's more fleshed out.

To install
```
make install
```

<h2>Document Indexer</h2>

```python
from suffix_array import SuffixArray

docs = [
    "The quick brown fox jumps over the lazy dog",
    "I am going to the store to buy some milk",
    "Uhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh"
]

## Creates and builds the index from documents.
suffix_array = SuffixArray(documents=docs, max_suffix_length=32)

## Query the index.
## Returns the documents that contain the query substring in a list.
records = suffix_array.query_records("the quick brown fox")
```

<h2>CSV Indexer</h2>
<p>Indexes and memory maps the file. Keeps only the suffix arrays for the text in memory (4 * N) where N = num_text_chars in search column.</p>

```python
from suffix_array import SuffixArray

## Creates and builds the index from csv
CSV_FILE      = "company_data.csv"
SEARCH_COLUMN = "company_name"

suffix_array = SuffixArray(
    csv_file=CSV_FILE,
    search_column=SEARCH_COLUMN,
    max_suffix_length=32
)

## Query the index.
## Returns the documents that contain the query substring in a list.
## (Returns all columns in list of dictionary (json records) format.)
records = suffix_array.query_records("netflix")
```

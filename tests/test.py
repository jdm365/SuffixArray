import polars as pl
from tqdm import tqdm
from time import perf_counter
import os

import pysubstringsearch
from suffix_array import SuffixArray


def test_pysubstring_search(csv_file: str, search_col: str):
    companies = pl.read_csv(
            csv_file,
            columns=[search_col],
            ).to_series(0).str.to_uppercase()

    N = 10_000
    rand_companies = [x for x in companies.sample(N)]

    init = perf_counter()
    writer = pysubstringsearch.Writer(
            index_file_path='tmp.idx'
            )

    for _, name in enumerate(tqdm(companies, desc='Building index')):
        writer.add_entry(name)

    writer.finalize()
    index_creation_time = perf_counter() - init

    reader = pysubstringsearch.Reader(
        index_file_path='tmp.idx'
        )

    times = []
    num_results = []

    init = perf_counter()
    for company in rand_companies:
        _init = perf_counter()
        vals = reader.search(company)
        time = (perf_counter() - _init) * 1e6
        times.append(time)
        num_results.append(len(vals))

    time = perf_counter() - init
    microseconds = 1e6 * time / 1000

    print("\n")
    print(f'Building pysubstringsearch index took  {index_creation_time:.2f} seconds')
    print(f'Querying pysubstringsearch took avg    {microseconds:.2f} microseconds')
    print(f'Querying pysubstringsearch took median {sorted(times)[len(times) // 2]:.2f} microseconds')
    print("\n\n\n\n")
    
    ## plot_query_times(times, num_results, outlier_stddev_threshold=3)


def test_suffix_array_document_constructor(csv_file: str, search_col: str):
    companies = pl.read_csv(
            csv_file,
            columns=[search_col],
            ).to_series(0).str.to_uppercase()

    N = 10_000
    rand_companies = [x for x in companies.sample(N)]
    list_companies = companies.to_list()

    init = perf_counter()
    suffix_array = SuffixArray(
            documents=list_companies, 
            max_suffix_length=32
            )
    index_creation_time = perf_counter() - init

    times = []
    num_results = []

    init = perf_counter()
    for company in rand_companies:
        _init = perf_counter()
        vals = suffix_array.query_records_2(company)
        print(vals)
        time = (perf_counter() - _init) * 1e6
        times.append(time)
        num_results.append(len(vals))

    time = perf_counter() - init
    microseconds = 1e6 * time / 1000
    
    print("\n")
    print(f'Building suffix array index took  {index_creation_time:.2f} seconds')
    print(f'Querying suffix array took avg    {microseconds:.2f} microseconds')
    print(f'Querying suffix array took median {sorted(times)[len(times) // 2]:.2f} microseconds')
    print("\n\n\n\n")

    ## plot_query_times(times, num_results, outlier_stddev_threshold=3)


def plot_query_times(query_times, number_of_results, outlier_stddev_threshold=3):
    ## Plot query time vs number of results
    import matplotlib.pyplot as plt
    import numpy as np

    ## Remove outlier number of results records
    number_of_results = np.array(number_of_results)
    max_num_results = np.mean(number_of_results) + outlier_stddev_threshold * np.std(number_of_results)
    mask = (number_of_results < max_num_results)
    query_times = np.array(query_times)[mask]
    number_of_results = number_of_results[mask]

    plt.scatter(query_times, number_of_results, alpha=0.5)

    ## Draw trendline
    z = np.polyfit(query_times, number_of_results, 1)
    p = np.poly1d(z)
    plt.plot(query_times, p(query_times), 'r--')

    ## Add grid
    plt.grid(True)

    ## Add legend
    plt.legend(['Data', 'Trendline'])

    plt.xlabel('Query time (microseconds)')
    plt.ylabel('Number of results')
    plt.title('Query time vs number of results')
    plt.show()



if __name__ == '__main__':
    DATA_DIR = '/home/jdm365/SearchApp/data'

    FILENAME = 'companies_sorted_100M.csv'
    ## FILENAME = 'companies_sorted.csv'
    ## FILENAME = 'companies_sorted_1M.csv'
    FILEPATH = os.path.join(DATA_DIR, FILENAME)

    test_pysubstring_search(FILEPATH, 'name')
    test_suffix_array_document_constructor(FILEPATH, 'name')

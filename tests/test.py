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
            ).to_series(0).fill_null('').str.to_uppercase()

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

    return times, num_results


def test_suffix_array_document_constructor(csv_file: str, search_col: str):
    companies = pl.read_csv(
            csv_file,
            columns=[search_col],
            ).to_series(0).fill_null('').str.to_uppercase()

    N = 10_000
    rand_companies = [x for x in companies.sample(N)]
    list_companies = companies.to_list()

    suffix_array = SuffixArray(
            max_suffix_length=32
            )
    init = perf_counter()
    suffix_array.construct_truncated_suffix_array_documents(documents=list_companies)
    index_creation_time = perf_counter() - init

    times = []
    num_results = []

    init = perf_counter()
    for company in rand_companies:
        _init = perf_counter()
        vals = suffix_array.query_records(company)
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

    plot_query_times(times, num_results, outlier_stddev_threshold=3)


def test_suffix_array_csv_constructor(csv_file: str, search_col: str):
    companies = pl.read_csv(
            csv_file,
            columns=[search_col],
            ).to_series(0).fill_null('').str.to_uppercase()

    N = 10_000
    rand_companies = [x for x in companies.sample(N)]

    suffix_array = SuffixArray(
            max_suffix_length=32
            )
    init = perf_counter()
    suffix_array.construct_truncated_suffix_array_from_csv(
            csv_file,
            search_col
            )
    index_creation_time = perf_counter() - init

    times = []
    num_results = []

    init = perf_counter()
    for company in tqdm(rand_companies):
        _init = perf_counter()
        vals = suffix_array.query_records(company)
        time = (perf_counter() - _init) * 1e6
        times.append(time)
        num_results.append(len(vals))

    time = perf_counter() - init
    microseconds = 1e6 * time / 1000
    print(sum(num_results) / len(num_results))
    
    print("\n")
    print(f'Building suffix array index took  {index_creation_time:.2f} seconds')
    print(f'Querying suffix array took avg    {microseconds:.2f} microseconds')
    print(f'Querying suffix array took median {sorted(times)[len(times) // 2]:.2f} microseconds')
    print("\n\n\n\n")

    ## plot_query_times(times, num_results, outlier_stddev_threshold=3)

    return times, num_results


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


def plot_query_times_multiple(
        query_times_sa, number_of_results_sa,
        query_times_pss, number_of_results_pss,
        outlier_stddev_threshold=3
        ):
    import matplotlib.pyplot as plt
    import numpy as np

    ## Remove outlier number of results records
    number_of_results_sa = np.array(number_of_results_sa)
    max_num_results = np.mean(number_of_results_sa) + outlier_stddev_threshold * np.std(number_of_results_sa)
    mask = (number_of_results_sa < max_num_results)
    query_times_sa = np.array(query_times_sa)[mask]
    number_of_results_sa = number_of_results_sa[mask]

    number_of_results_pss = np.array(number_of_results_pss)
    max_num_results = np.mean(number_of_results_pss) + outlier_stddev_threshold * np.std(number_of_results_pss)
    mask = (number_of_results_pss < max_num_results)
    query_times_pss = np.array(query_times_pss)[mask]
    number_of_results_pss = number_of_results_pss[mask]

    plt.scatter(query_times_sa, number_of_results_sa, alpha=0.5)
    plt.scatter(query_times_pss, number_of_results_pss, alpha=0.5)

    ## Draw trendlines
    z = np.polyfit(query_times_sa, number_of_results_sa, 1)
    p = np.poly1d(z)
    plt.plot(query_times_sa, p(query_times_sa), 'r--')

    z = np.polyfit(query_times_pss, number_of_results_pss, 1)
    p = np.poly1d(z)
    plt.plot(query_times_pss, p(query_times_pss), 'g--')

    ## Add grid
    plt.grid(True)

    ## Add legend
    ## plt.legend(['Data (SA)', 'Trendline (SA)', 'Data (PSS)', 'Trendline (PSS)'])
    plt.legend([
        'Suffix Array',
        'PySubstringSearch',
        'Suffix Array Trendline',
        'PySubstringSearch Trendline'
        ])

    plt.xlabel('Query time (microseconds)')
    plt.ylabel('Number of results')
    plt.title('Query time vs number of results')
    plt.show()


if __name__ == '__main__':
    ## DATA_DIR = '/home/jdm365/SearchApp/data'
    CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
    DATA_DIR = '../data'
    DATA_DIR = os.path.join(CURRENT_DIR, DATA_DIR)

    ## FILENAME = 'companies_sorted_100M.csv'
    ## FILENAME = 'companies_sorted.csv'
    ## FILENAME = 'companies_sorted_1M.csv'
    FILENAME = 'companies-2023-q4-sm.csv'
    FILEPATH = os.path.join(DATA_DIR, FILENAME)

    ## test_suffix_array_document_constructor(FILEPATH, 'name')
    times_sa, num_results_sa = test_suffix_array_csv_constructor(FILEPATH, 'name')
    times_pss, num_results_pss = test_pysubstring_search(FILEPATH, 'name')
    plot_query_times_multiple(times_sa, num_results_sa, times_pss, num_results_pss)

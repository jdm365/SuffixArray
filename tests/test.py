import pandas as pd
from tqdm import tqdm
from time import perf_counter
import os

import pysubstringsearch
from suffix_array import SuffixArray



if __name__ == '__main__':
    DATA_DIR = '/home/jdm365/SearchApp/data'
    FILENAME = 'companies_sorted_100M.csv'
    ## FILENAME = 'companies_sorted.csv'
    ## FILENAME = 'companies_sorted_1M.csv'
    FILEPATH = os.path.join(DATA_DIR, FILENAME)

    companies = pd.read_csv(
            FILEPATH,
            usecols=['name'],
            ## nrows=1000
            engine='pyarrow'
            )['name'].str.upper()
    '''
    with open(os.path.join(DATA_DIR, 'names_only.txt'), 'r') as f:
        companies = f.readlines()
    '''

    N = 10_000
    ## rand_companies = companies.sample(1000)
    rand_companies = [x[:8] for x in companies.sample(N)]

    init = perf_counter()
    writer = pysubstringsearch.Writer(
            index_file_path='tmp.idx'
            )

    for idx, name in enumerate(tqdm(companies, desc='Building index')):
        writer.add_entry(name)
    ## joblib.dump(self.inverted_index, 'inverted_index.joblib', protocol=4)

    writer.finalize()
    print(f'Building pysubstringsearch index took {perf_counter() - init:.2f} seconds')

    reader = pysubstringsearch.Reader(
        index_file_path='tmp.idx'
        )

    init = perf_counter()
    idxs = reader.search('NETFLIX')
    print(f'Querying pysubstringsearch index took {1e6 * (perf_counter() - init):.2f} microseconds')

    init = perf_counter()
    for company in rand_companies:
        idxs = reader.search(company)
    time = perf_counter() - init
    microseconds = 1e6 * time / 1000
    print(f'Querying pysubstringsearch index took avg {microseconds:.2f} microseconds')


    init = perf_counter()
    suffix_array = SuffixArray(documents=companies, max_suffix_length=32)
    print(f'Building suffix array index took {perf_counter() - init:.2f} seconds')

    init = perf_counter()
    idxs = suffix_array.query('NETFLIX')
    print(f'Querying suffix array took {1e6 * (perf_counter() - init):.2f} microseconds')

    print(companies.iloc[idxs])

    times = []
    num_results = []

    init = perf_counter()
    for company in tqdm(rand_companies, desc='Querying suffix array'):
        _init = perf_counter()
        idxs = suffix_array.query(company)
        time = (perf_counter() - _init) * 1e6
        ## print(f'Querying {company} took {time:.2f} microseconds\n')
        times.append(time)
        num_results.append(len(idxs))

    time = perf_counter() - init
    microseconds = 1e6 * time / 1000
    print(f'Querying suffix array took avg    {microseconds:.2f} microseconds')
    print(f'Querying suffix array took median {sorted(times)[len(times) // 2]:.2f} microseconds')

    os.system(f'du -sh {FILEPATH}')
    os.system('du -sh suffix_array_data')
    os.system('du -sh tmp.idx')

    os.remove('tmp.idx')
    ## os.system('rm -rf suffix_array_data')


    '''
    ## Plot query time vs number of results
    import matplotlib.pyplot as plt
    import numpy as np

    plt.scatter(times, num_results, alpha=0.5)

    ## Draw trendline
    z = np.polyfit(times, num_results, 1)
    p = np.poly1d(z)
    plt.plot(times, p(times), 'r--')

    ## Add grid
    plt.grid(True)

    ## Add legend
    plt.legend(['Data', 'Trendline'])

    plt.xlabel('Query time (microseconds)')
    plt.ylabel('Number of results')
    plt.title('Query time vs number of results')
    plt.show()
    '''

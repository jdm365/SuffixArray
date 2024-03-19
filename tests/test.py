import pandas as pd
from tqdm import tqdm
from time import perf_counter
import os

import pysubstringsearch
from suffix_array import SuffixArray



if __name__ == '__main__':
    DATA_DIR = '/home/jdm365/SearchApp/data'
    FILENAME = 'companies_sorted.csv'
    FILEPATH = os.path.join(DATA_DIR, FILENAME)

    companies = pd.read_csv(
            FILEPATH,
            usecols=['name'],
            ## nrows=1000
            )['name'].str.upper()
    '''
    with open(os.path.join(DATA_DIR, 'names_only.txt'), 'r') as f:
        companies = f.readlines()
    '''

    rand_companies = companies.sample(1000)

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

    os.remove('tmp.idx')

    init = perf_counter()
    for company in rand_companies:
        idxs = reader.search(company)
    time = perf_counter() - init
    microseconds = 1e6 * time / 1000
    print(f'Querying pysubstringsearch index took avg {microseconds:.2f} microseconds')


    init = perf_counter()
    suffix_array = SuffixArray(documents=companies, max_suffix_length=64)
    print(f'Building suffix array index took {perf_counter() - init:.2f} seconds')

    init = perf_counter()
    idxs = suffix_array.query('NETFLIX')
    print(f'Querying suffix array took {1e6 * (perf_counter() - init):.2f} microseconds')

    print(companies.iloc[idxs])

    init = perf_counter()
    for company in tqdm(rand_companies, desc='Querying suffix array'):
        idxs = suffix_array.query(company)
    time = perf_counter() - init
    microseconds = 1e6 * time / 1000
    print(f'Querying suffix array took avg {microseconds:.2f} microseconds')

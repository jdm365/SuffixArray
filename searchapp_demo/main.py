from flask import Flask, request, jsonify
from flask_cors import CORS

import os
import polars as pl
from time import perf_counter
from suffix_array import SuffixArrayEngine
import json
import csv
from functools import lru_cache


app = Flask(__name__)
CORS(app)

@app.route('/search', methods=['GET'])
def search():
    query = request.args.get('query', '')

    init = perf_counter()

    results = search_app.get_search_results(query)

    time_taken_ms = int(1e3) * (perf_counter() - init)
    return jsonify({'results': results, 'time_taken_ms': time_taken_ms})


class SearchApp:
    def __init__(self, csv_filename: str) -> None:

        self.search_col = 'name'

        self.csv_filename = csv_filename
        self.column_names = self.get_column_names()

        self.reader = SuffixArrayEngine(max_suffix_length=32)

        '''
        init = perf_counter()
        self.reader.construct_truncated_suffix_array_from_csv(
                filename=self.csv_filename,
                search_column="name"
                )
        print(f'Building suffix array took {perf_counter() - init:.2f} seconds')
        '''
        init = perf_counter()
        self.reader.construct_truncated_suffix_array_from_csv(
                filename=self.csv_filename,
                search_column=self.search_col
                )
        print(f'Building suffix array took {perf_counter() - init:.2f} seconds')

        '''
        init = perf_counter()
        self.reader.save(save_dir='suffix_array_data')
        print(f'Saving suffix array took {perf_counter() - init:.2f} seconds')
        sys.stdout.flush()

        del self.reader
        
        self.reader = SuffixArray(max_suffix_length=32)

        init = perf_counter()
        self.reader.load(save_dir='suffix_array_data')
        print(f'Loading suffix array took {perf_counter() - init:.2f} seconds')
        sys.stdout.flush()
        '''


    def get_column_names(self):
        if self.csv_filename.endswith('.csv'):
            with open(self.csv_filename, newline='') as csvfile:
                reader = csv.reader(csvfile)
                columns = next(reader)

        elif self.csv_filename.endswith('.json'):

            ## Read first line
            with open(self.csv_filename, 'r') as f:
                data = f.readlines()
                columns = list(json.loads(data[0]).keys())
        else:
            raise ValueError(f"Unknown file type: {self.csv_filename}\n Please provide a .csv or .json file.")

        cols = [x.lower() for x in columns if x.strip() != '']
        return cols
            

    @lru_cache(maxsize=16)
    def get_search_results(self, query: str) -> list:
        query = query.lower()
        if len(query) == 0:
            return []
        init = perf_counter()
        print(f'Query: {query}')

        init = perf_counter()
        vals = self.reader.query_records(query, k=1000)

        print(f'Query took {1000 * (perf_counter() - init):.2f} milliseconds')

        if len(vals) == 0:
            return []

        return vals




if __name__ == '__main__':
    ## DATA_DIR = '/home/jdm365/SearchApp/data'
    CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
    DATA_DIR = CURRENT_DIR
    ## FILENAME = 'companies_sorted_1M.csv'
    ## FILENAME = 'people_data_labs_sample.csv'
    FILENAME = '../data/companies-2023-q4-sm.csv'
    ## FILENAME = 'companies_100M.csv'
    FILEPATH = os.path.join(DATA_DIR, FILENAME)

    search_app = SearchApp(
            csv_filename=FILEPATH,
            )

    os.system(f"open {os.path.join(CURRENT_DIR, 'index.html')}")

    app.run()

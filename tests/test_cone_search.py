import unittest
from unittest.mock import patch, MagicMock
import pandas as pd
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
import sys
import os

# Add the module path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../mindless_scripts/astro_tools')))

from cone_search import perform_radio_search, perform_simbad_search

class TestConeSearch(unittest.TestCase):

    def setUp(self):
        self.coord = SkyCoord(ra=10.0, dec=20.0, unit='deg', frame='icrs')
        self.radius = 10.0 # arcsec

    @patch('cone_search.Vizier')
    def test_perform_radio_search(self, mock_vizier):
        # Mock Vizier response
        mock_v = MagicMock()
        mock_vizier.return_value = mock_v
        
        # Create a dummy result table
        data = {
            '_RAJ2000': [10.0001],
            '_DEJ2000': [20.0001],
            '_r': [0.5] # arcsec
        }
        mock_table = MagicMock()
        mock_table.to_pandas.return_value = pd.DataFrame(data)
        
        # Vizier returns a TableList, which is list-like
        mock_v.query_region.return_value = [mock_table]
        
        results = perform_radio_search(self.coord, self.radius)
        
        self.assertIsInstance(results, pd.DataFrame)
        self.assertFalse(results.empty)
        self.assertIn('ID', results.columns)
        self.assertIn('Catalogue', results.columns)
        # We expect at least one result per catalogue if we mocked it to return something for all
        # But perform_radio_search iterates catalogues. 
        # Our mock returns the same thing for every call.
        # So we should get len(RADIO_CATALOGUES) results.
        self.assertGreater(len(results), 0)

    @patch('cone_search.Simbad')
    def test_perform_simbad_search(self, mock_simbad):
        # Mock Simbad response
        mock_s = MagicMock()
        mock_simbad.return_value = mock_s
        
        data = {
            'MAIN_ID': ['Test Source'],
            'RA_d': [10.0002],
            'DEC_d': [20.0002],
            'PMRA': [10.0],
            'PMDEC': [-5.0],
            'PLX_VALUE': [2.0]
        }
        mock_table = MagicMock()
        mock_table.to_pandas.return_value = pd.DataFrame(data)
        
        mock_s.query_region.return_value = mock_table
        
        results = perform_simbad_search(self.coord, self.radius)
        
        self.assertIsInstance(results, pd.DataFrame)
        self.assertFalse(results.empty)
        self.assertEqual(results.iloc[0]['ID'], 'Test Source')
        self.assertIn('PMRA', results.columns)
        self.assertIn('PMDec', results.columns)
        self.assertAlmostEqual(results.iloc[0]['Distance_arcsec'], self.coord.separation(SkyCoord(10.0002, 20.0002, unit='deg')).arcsec, places=1)

    @patch('cone_search.Vizier')
    def test_perform_radio_search_no_results(self, mock_vizier):
        mock_v = MagicMock()
        mock_vizier.return_value = mock_v
        mock_v.query_region.return_value = [] # No results
        
        results = perform_radio_search(self.coord, self.radius)
        self.assertTrue(results.empty)

    @patch('cone_search.Simbad')
    def test_perform_simbad_search_no_results(self, mock_simbad):
        mock_s = MagicMock()
        mock_simbad.return_value = mock_s
        mock_s.query_region.return_value = None
        
        results = perform_simbad_search(self.coord, self.radius)
        self.assertTrue(results.empty)

if __name__ == '__main__':
    unittest.main()

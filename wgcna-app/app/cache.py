# cache.py

import os
import glob
import scanpy as sc
import threading
from config import Config

class AnnDataCache:
    """
    Singleton class to manage caching of AnnData objects loaded from .h5ad files.
    """

    _instance = None
    _lock = threading.Lock()

    def __new__(cls, h5ad_dir):
        """
        Ensure only one instance of AnnDataCache exists (Singleton pattern).
        """
        
        if cls._instance is None:
            with cls._lock:
                if cls._instance is None:
                    cls._instance = super(AnnDataCache, cls).__new__(cls)
                    cls._instance.h5ad_dir = h5ad_dir
                    cls._instance.cache = {}
                    cls._instance._load_cache()
        return cls._instance

    def _load_cache(self):
        """
        Load specific .h5ad files based on APP_NAME into the cache.
        """

        for adata_file in glob.glob(os.path.join(self.h5ad_dir, "*.h5ad")):
            key = os.path.basename(adata_file).split('.')[0]

            if Config.APP_NAME == "Ceratopteris" and key != "CR":
                continue  # Skip all except "CR"
            elif Config.APP_NAME == "RanOmics" and key == "CR":
                continue  # Skip "CR" for RanOmics

            print(f"Loading data from {adata_file}")
            adata = sc.read_h5ad(adata_file)
            self.cache[key] = adata

    def get_adata(self, key):
        """
        Retrieve an AnnData object from the cache by key.

        Parameters:
        - key (str): The key corresponding to the desired AnnData object.

        Returns:
        - AnnData or None: The cached AnnData object or None if not found.
        """

        return self.cache.get(key)

    def keys(self):
        """
        Get all keys in the cache.

        Returns:
        - list: List of keys.
        """

        return list(self.cache.keys())


def get_adata_cache():
    """"
    Get the singleton instance of AnnDataCache.

    Returns:
    - AnnDataCache: The AnnDataCache instance.
    """
    
    if AnnDataCache._instance is None:
        AnnDataCache(Config.H5AD_DIR)
        
    return AnnDataCache._instance

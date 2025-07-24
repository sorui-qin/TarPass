'''
Author: Rui Qin
Date: 2025-06-24 14:50:13
LastEditTime: 2025-07-23 14:17:57
Description: 
'''
from rdkit.Chem import Mol, FilterCatalog
from rdkit.Chem.FilterCatalog import FilterCatalogParams

def catalog(catalog_type):
    params = FilterCatalogParams()
    params.AddCatalog(catalog_type)
    return FilterCatalog.FilterCatalog(params)

_PAINS = catalog(FilterCatalogParams.FilterCatalogs.PAINS)
_SURECHEMBL = catalog(FilterCatalogParams.FilterCatalogs.CHEMBL_SureChEMBL)
_GLAXO = catalog(FilterCatalogParams.FilterCatalogs.CHEMBL_Glaxo)

def match_pains(mol:Mol):
    """Pan Assay Interference Compounds (PAINS) Filters.  
    Ref: `J. Med. Chem. 2010, 53, 7, 2719–2740`
    """
    return int(_PAINS.HasMatch(mol))

def match_SureChEMBL(mol:Mol):
    """Structural alerts such toxicophores appiled in SureChEMBL, based on ToxAlerts.
    Ref: `J. Chem. Inf. Model. 2012, 52, 8, 2310–2316`
    """
    return int(_SURECHEMBL.HasMatch(mol))

def match_glaxo(mol:Mol):
    """Hard Filters raised by Glaxo Wellcome.  
    Containing reactive functional groups, unsuitable leads, and unsuitable natural products.  
    Ref: `J. Chem. Inf. Comput. Sci. 1999, 39, 5, 897–902`
    """
    return int(_GLAXO.HasMatch(mol))
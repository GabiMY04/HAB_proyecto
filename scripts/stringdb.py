# stringdb.py
# Wrapper mínimo para consultar enriquecimiento en STRING y devolver un DataFrame.

from __future__ import annotations
import requests
import pandas as pd
from urllib.parse import quote_plus

BASE_URL = "https://string-db.org/api/json/enrichment"

def _normalize_genes(genes) -> list[str]:
    if genes is None:
        return []
    if isinstance(genes, str):
        genes = [g.strip() for chunk in genes.split(",") for g in chunk.splitlines()]
    return [g for g in genes if g]

def get_enrichment(genes, species: int = 3702) -> pd.DataFrame:
    ids = _normalize_genes(genes)
    if not ids:
        return pd.DataFrame()
    identifiers = "%0d".join(quote_plus(g) for g in ids)  # separador que usa STRING
    r = requests.get(BASE_URL, params={"identifiers": identifiers, "species": species}, timeout=60)
    r.raise_for_status()
    df = pd.DataFrame(r.json())
    if df.empty:
        return df
    if "term" not in df.columns and "term_id" in df.columns:
        df["term"] = df["term_id"]
    if "description" not in df.columns and "term_desc" in df.columns:
        df["description"] = df["term_desc"]
    if "fdr" not in df.columns:
        for cand in ("FDR", "fdr_value", "adj_p", "Adjusted P-value"):
            if cand in df.columns:
                df["fdr"] = df[cand]; break
        else:
            df["fdr"] = pd.NA
    if "inputGenes" not in df.columns:
        df["inputGenes"] = [[] for _ in range(len(df))]
    try:
        df = df.sort_values("fdr", ascending=True)
    except Exception:
        pass
    return df

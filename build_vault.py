#!/usr/bin/env python3
"""Parse references.bib and build Obsidian vault notes with PMID and abstract."""

import re
import time
import xml.etree.ElementTree as ET
from urllib.parse import quote_plus
from urllib.request import urlopen
import json

BIB_FILE = "lit/references.bib"
VAULT_DIR = "vault"

def parse_bib(bib_path):
    """Yield dicts for each @article entry."""
    entries = []
    current = None
    field_re = re.compile(r'^\s+(\w+)\s*=\s*\{(.+?)\}\s*(?:,|$)')
    key_re = re.compile(r'^@\w+\{(\w+),')
    with open(bib_path) as f:
        for line in f:
            m = key_re.match(line)
            if m:
                if current:
                    entries.append(current)
                current = {"key": m.group(1)}
                continue
            if current is None:
                continue
            m = field_re.match(line)
            if m:
                current[m.group(1).lower()] = m.group(2)
            if line.strip() == '}':
                entries.append(current)
                current = None
    if current:
        entries.append(current)
    return entries

def fetch_pubmed_data(doi):
    """Return (pmid, abstract) for a DOI via NCBI E-utilities."""
    # Step 1: search for PMID by DOI
    search_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        f"?db=pubmed&term={quote_plus(doi)}[DOI]&retmode=json"
    )
    try:
        with urlopen(search_url, timeout=15) as resp:
            search_data = json.loads(resp.read())
    except Exception:
        return None, None
    id_list = search_data.get("esearchresult", {}).get("idlist", [])
    if not id_list:
        return None, None
    pmid = id_list[0]
    # Step 2: fetch abstract
    fetch_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
        f"?db=pubmed&id={pmid}&retmode=xml&rettype=abstract"
    )
    try:
        with urlopen(fetch_url, timeout=15) as resp:
            xml_data = resp.read()
    except Exception:
        return pmid, None
    root = ET.fromstring(xml_data)
    abstract = ""
    for abstract_elem in root.iter("AbstractText"):
        label = abstract_elem.get("Label", "")
        text = "".join(abstract_elem.itertext())
        if label:
            abstract += f"**{label}:** {text}\n\n"
        else:
            abstract += text + "\n\n"
    abstract = abstract.strip()
    return pmid, abstract if abstract else None

def make_filename(key):
    """Convert citation key to a clean markdown filename."""
    return re.sub(r'[<>:"/\\|?*]', '', key).strip() + ".md"

def clean_abstract_for_obsidian(text):
    """Remove problematic characters for Obsidian."""
    if text is None:
        return ""
    text = re.sub(r'[\x00-\x08\x0b\x0c\x0e-\x1f]', '', text)
    return text

def main():
    entries = parse_bib(BIB_FILE)
    print(f"Found {len(entries)} entries in .bib file")

    for i, entry in enumerate(entries):
        key = entry.get("key", f"entry_{i}")
        title = entry.get("title", "").replace("{", "").replace("}", "")
        doi = entry.get("doi", "")
        authors = entry.get("author", "")
        journal = entry.get("journal", "")
        year = entry.get("year", "")
        note = entry.get("note", "")

        filename = make_filename(key)
        filepath = f"{VAULT_DIR}/{filename}"

        print(f"\n[{i+1}/{len(entries)}] {key}: fetching PMID/abstract...")
        pmid, abstract = None, None
        if doi:
            pmid, abstract = fetch_pubmed_data(doi)
            time.sleep(0.4)  # NCBI rate limit
        else:
            print("  No DOI, skipping PubMed lookup")

        lines = []
        lines.append(f"# {title}")
        lines.append("")
        lines.append(f"**Citation key:** `{key}`")
        if authors:
            lines.append(f"**Authors:** {authors}")
        if journal:
            lines.append(f"**Journal:** {journal}")
        if year:
            lines.append(f"**Year:** {year}")
        if doi:
            lines.append(f"**DOI:** [{doi}](https://doi.org/{doi})")
        if pmid:
            lines.append(f"**PMID:** [{pmid}](https://pubmed.ncbi.nlm.nih.gov/{pmid}/)")
        else:
            lines.append("**PMID:** Not found in PubMed")
        if note:
            lines.append(f"**Note:** {note}")
        lines.append("")
        lines.append("---")
        lines.append("")

        if abstract:
            lines.append("## Abstract")
            lines.append("")
            lines.append(clean_abstract_for_obsidian(abstract))
        else:
            lines.append("## Abstract")
            lines.append("")
            lines.append("*No abstract available from PubMed.*")

        if doi and not pmid:
            lines.append("")
            lines.append(f"> **PubMed lookup note:** The DOI `{doi}` did not return a PubMed record. "
                         "The article may not be indexed in PubMed.")

        content = "\n".join(lines)
        with open(filepath, "w") as f:
            f.write(content)
        print(f"  Wrote {filepath} (PMID: {pmid})")

if __name__ == "__main__":
    main()

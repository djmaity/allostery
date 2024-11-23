import glob
import os

import pandas as pd
import tqdm
import xmltodict


def parse_allosteric_site(alloteric_site):
    """Convert the list of residues in ASD from "Chain A:HIS25,TYR258; Chain B:VAL325" to ['A-HIS-25', 'A-TYR-258', 'B-VAL-325']"""
    residues = []
    for chain_string in alloteric_site.split('; '):
        chain_name, residue_string = chain_string.split(':')
        chain_id = chain_name[-1]
        for residue in residue_string.split(','):
            res_name, res_id = residue[:3], residue[3:]
            residues.append(f'{chain_id}-{res_name}-{res_id}')
    return residues


def parse_asd_xml(asd_xml_dir):
    data = []
    for xml_file in tqdm.tqdm(glob.glob(os.path.join(asd_xml_dir, '*.xml'))):
        xml_string = open(xml_file, 'r').read()
        xml_string = xml_string.replace('&#x2;', '')  # Remove invalid XML character
        protein = xmltodict.parse(xml_string)
    
        target_gene = ''
        target_id = protein['Organism_Record']['Organism_ID']
    
        if 'Gene' in protein['Organism_Record'] and 'Gene_Name' in protein['Organism_Record']['Gene']:
            target_gene = protein['Organism_Record']['Gene']['Gene_Name']
    
        organism = protein['Organism_Record']['Organism']
    
        # Check if the list of allosteric sites in present in the XML file
        if 'Allosteric_Site_List' in protein['Organism_Record']:
            allosteric_sites = protein['Organism_Record']['Allosteric_Site_List']['Allosteric_Site']
    
            # Check if more than one allosteric site is present
            if isinstance(allosteric_sites, list):
                allosteric_sites = allosteric_sites
            else:
                allosteric_sites = [allosteric_sites]
        else:
            allosteric_sites = []
    
        for site in allosteric_sites:
            if site:
                pdb_uniprot = site['PDB_UniProt_ID'] if 'PDB_UniProt_ID' in site else ''
                allosteric_pdb = site['Allosteric_PDB']
                modulator_serial = site['Modulator_ASD_ID']
                modulator_alias = site['Modulator_Alias']
                modulator_chain = site['Modulator_Chain']
                modulator_class = site['Modulator_Class'] if 'Modulator_Class' in site else ''
                modulator_feature = site['Modulator_Feature'] if 'Modulator_Feature' in site else ''
                modulator_name = site['Modulator_Name'] if 'Modulator_Name' in site else ''
                modulator_resi = site['Modulator_Residue'] if 'Modulator_Residue' in site else ''
                function = site['Function'] if 'Function' in site else ''
                position = site['Position'] if 'Position' in site else ''
                pubmed_id = site['PubMed_ID'] if 'PubMed_ID' in site else ''
                ref_title = site['PubMed_Title'] if 'PubMed_Title' in site else ''
                site_overlap = site['Site_Overlap'] if 'Site_Overlap' in site else ''
                allosteric_site_residue = parse_allosteric_site(site['Allosteric_Site_Residue']) if 'Allosteric_Site_Residue' in site else []
    
                data.append([
                    target_id, target_gene, organism, pdb_uniprot, allosteric_pdb,
                    modulator_serial, modulator_alias, modulator_chain,
                    modulator_class, modulator_feature, modulator_name,
                    modulator_resi, function, position, pubmed_id, ref_title,
                    site_overlap, allosteric_site_residue])
    
    return pd.DataFrame(data, columns=[
        'target_id', 'target_gene', 'organism', 'pdb_uniprot', 'allosteric_pdb',
        'modulator_serial', 'modulator_alias', 'modulator_chain', 'modulator_class',
        'modulator_feature', 'modulator_name', 'modulator_resi', 'function',
        'position', 'pubmed_id', 'ref_title', 'site_overlap',
        'allosteric_site_residue'])


if __name__ == "__main__":
    # Extract the ASD_Release_202306_XF.tar.gz archive obtained from
    # https://mdl.shsmu.edu.cn/ASD/module/download/download.jsp?tabIndex=1
    # to the ../data directory
    df_asd = parse_asd_xml('../data/ASD_Release_202306_XF/')

    # Write the CSV file from the parsed data
    df_asd.to_csv('../output/ASD_Release_202306.csv', index=False)
import pandas as pd
import requests
import os
import tqdm
import urllib
import Bio.PDB


def get_pdb_data(pdb_ids):
    query = """query structure($pdb_ids: [String!]!) {
  entries(entry_ids: $pdb_ids) {
    rcsb_id
    assemblies {
      rcsb_id
      polymer_entity_instances {
        rcsb_id
        polymer_entity {
          entity_poly {
            pdbx_seq_one_letter_code_can
          }
          uniprots {
            rcsb_id
          }
        }
      }
      rcsb_struct_symmetry {
        kind
        oligomeric_state
        stoichiometry
      }
    }
  }
}"""
    response = requests.post(url='https://data.rcsb.org/graphql', json={"query": query, "variables": {"pdb_ids": list(pdb_ids)}})
    response_data = response.json()

    pdb_data = []
    for entry in response_data['data']['entries']:
        pdb_id = entry['rcsb_id']

        first_assembly = entry['assemblies'][0]
        # assembly_id = first_assembly['rcsb_id'].split('-')[1]

        chain_uniprot_mapping = []
        for polymer_entity in first_assembly['polymer_entity_instances']:
            chain = polymer_entity['rcsb_id'].split('.')[1]
            sequence = polymer_entity['polymer_entity']['entity_poly']['pdbx_seq_one_letter_code_can']
            polymer_entity_uniprot = polymer_entity['polymer_entity']['uniprots']
            uniprot_ids = []
            if polymer_entity_uniprot:
                for uniprot_id in polymer_entity_uniprot:
                    uniprot_ids.append(uniprot_id['rcsb_id'])
            chain_uniprot_mapping.append([chain, uniprot_ids, sequence])

        if first_assembly['rcsb_struct_symmetry'] is not None:
            for symmetry in first_assembly['rcsb_struct_symmetry']:
                if symmetry['kind'] == 'Global Symmetry':
                    oligomeric_state = symmetry['oligomeric_state']
                    stoichiometry = symmetry['stoichiometry']
        else:
            oligomeric_state = None
            stoichiometry = None
        pdb_data.append([pdb_id, chain_uniprot_mapping, oligomeric_state, stoichiometry])

    return pd.DataFrame(pdb_data, columns=['pdb_id', 'chain_uniprot_mapping', 'oligomeric_state', 'stoichiometry'])

####################################################################################

def uniprot_from_pdb(pdb_ids):
    """ Get the UniProt IDs for the chains in the structures with the PDB IDs """
    pdb_ids_string = '", "'.join(pdb_ids)

    body = 'query {entries(entry_ids: ["' + pdb_ids_string + '"])' + """
      {
        rcsb_id
        polymer_entities {
          uniprots {
            rcsb_id
          }
        }
      }
    }
    """
    response = requests.post(url='https://data.rcsb.org/graphql', json={"query": body})
    response_data = response.json()

    id_mapping = []
    for record in response_data['data']['entries']:
        rcsb_id = record['rcsb_id']
        uniprot_ids = set()
        for entity in record['polymer_entities']:
            if entity['uniprots']:
                uniprot_ids.add(entity['uniprots'][0]['rcsb_id'])
        id_mapping.append([rcsb_id, uniprot_ids])
        
    return pd.DataFrame(id_mapping, columns=['pdb_id', 'uniprot_id_from_pdb'])


def uniprot_pdb_fix(df):
    df_copy = df[['pdb_id', 'uniprot_id']].drop_duplicates().copy()

    # Get the UniProt IDs of the chains in the proteins from PDB website
    df_pdb = uniprot_from_pdb(df_copy['pdb_id'])
    
    # Count the number of UniProt IDs obtained from PDB website for each PDB ID
    df_pdb['num_uniprot'] = df_pdb['uniprot_id_from_pdb'].apply(len)
    
    df_merged = df_copy.merge(df_pdb, how='left')

    # Check if the UniProt ID in ASD is present in the UniProt IDs obtained from PDB website
    df_merged['present'] = df_merged.apply(lambda x: x['uniprot_id'] in x['uniprot_id_from_pdb'], axis=1)
    df_present = df_merged.loc[df_merged['present'], ['pdb_id', 'uniprot_id']]

    # UniProt IDs absent in the ASD but has only one UniProt ID downloaded from PDB
    df_absent_single = df_merged[~df_merged['present'] & (df_merged['num_uniprot'] == 1)].copy()
    
    # print("PDB IDs Absent from ASD with Single UniProt ID fetched from the PDB")
    # display(df_absent_single)

    # Place the only one UniProt ID from the set in Uniprot ID column to pdb_uniprot column
    if not df_absent_single.empty:
        df_absent_single['uniprot_id'] = df_absent_single.apply(lambda x: next(iter(x['uniprot_id_from_pdb'])), axis=1)
    df_absent_single = df_absent_single[['pdb_id', 'uniprot_id']]
    
    df_absent_multiple = df_merged[~df_merged['present'] & (df_merged['num_uniprot'] != 1)]
    print("PDB IDs Absent from ASD with Zero or Multiple UniProt ID fetched from the PDB")
    display(df_absent_multiple[['pdb_id', 'uniprot_id', 'uniprot_id_from_pdb']].sort_values('pdb_id'))

    df_output = pd.concat([df_present, df_absent_single])
    return df_output.sort_values(['pdb_id', 'uniprot_id']).reset_index(drop=True)


# def replace_obsolete_pdbs(pdb_set, obsolete_pdbs):
#     for pdb in pdb_set & set(obsolete_pdbs.keys()):
#         pdb_set.remove(pdb)
#         pdb_set.add(obsolete_pdbs[pdb])
#         print(f'PDB ID {pdb} has been superseded by {obsolete_pdbs[pdb]}')




def download_structures(pdb_ids: list[str], download_dir: os.path = './') -> None:

    if not os.path.exists(download_dir):
        os.makedirs(download_dir)

    for pdb_id in tqdm.tqdm(pdb_ids):
        pdb_id = pdb_id.lower()
        pdb_file = os.path.join(download_dir, f'{pdb_id}.pdb')
        cif_file = os.path.join(download_dir, f'{pdb_id}.cif')
        if not os.path.exists(pdb_file) and not os.path.exists(cif_file):
            try:
                urllib.request.urlretrieve(f'https://files.rcsb.org/download/{pdb_id}.pdb1', pdb_file)
            except urllib.error.HTTPError:
                try:
                    urllib.request.urlretrieve(f'https://files.rcsb.org/download/{pdb_id}-assembly1.cif', cif_file)
                except urllib.error.HTTPError:
                    print('Could not download', pdb_id.upper())


#########################################################


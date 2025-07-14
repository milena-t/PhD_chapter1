import requests
import xml.etree.ElementTree as ET
from ete3 import NCBITaxa, Tree

accession_nos = ["GCA_000355655.1","GCA_000390285.2","GCA_000500325.2","GCA_000648695.2","GCA_000699045.2","GCA_001412225.1","GCA_002938485.2","GCA_003013835.2","GCA_004193795.1","GCA_008802855.1","GCA_011009095.1","GCA_014462685.1","GCA_015345945.1","GCA_020466585.2","GCA_020466635.2","GCA_027724725.1","GCA_027725215.1","GCA_029955315.1","GCA_963669975.1"]


# get taxonomic IDs and scientific names of all the assemblies from ncbi:
# actually running this takes 10 seconds or so, so I will not do that every time
if False:
    def search_assembly_id(assembly_accession):
        base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi'
        params = {
            'db': 'assembly',  # Use 'assembly' database
            'term': assembly_accession,
            'retmode': 'json'
        }
        try:
            response = requests.get(base_url, params=params)
            response.raise_for_status()  # Raise an exception for 4xx or 5xx status codes
            data = response.json()
            return data.get('esearchresult', {}).get('idlist', [])
        except requests.exceptions.RequestException as e:
            print("Error making request:", e)
            return []

    def get_tax_id(assembly_id):
        base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi'
        params = {
            'db': 'assembly',
            'id': assembly_id,
            'retmode': 'json'
        }
        try:
            response = requests.get(base_url, params=params)
            response.raise_for_status()  # Raise an exception for 4xx or 5xx status codes
            data = response.json()
            if 'result' in data and assembly_id in data['result']:
                assembly_info = data['result'][assembly_id]
                taxonomic_id = assembly_info.get('speciestaxid')
                species_name = assembly_info.get("speciesname")
                return [taxonomic_id, species_name]
            else:
                return None
        except requests.exceptions.RequestException as e:
            print("Error making request:", e)
            return None
        except ValueError as e:
            print("Error decoding JSON:", e)
            return None



    assembly_IDs = [search_assembly_id(ac_no)[0] for ac_no in accession_nos]
    IDs = [get_tax_id(as_id) for as_id in assembly_IDs]
    taxonomic_IDs = [id[0] for id in IDs]
    species_names = [id[1] for id in IDs]

# results from the above functions copy-pasted from the command line
taxonomic_IDs = ['77166', '217634', '7539', '166361', '224129', '110193', '7048', '50389', '1661398', '7054', '2038154', '354439', '41895', '77166', '77166', '2755281', '7067', '197179', '200917']
#species_names = ['Bruchidius siliquastri', 'Dendroctonus ponderosae', 'Anoplophora glabripennis', 'Leptinotarsa decemlineata', 'Onthophagus taurus', 'Agrilus planipennis', 'Nicrophorus vespilloides', 'Sitophilus oryzae', 'Diabrotica virgifera', 'Asbolus verrucosus', 'Photinus pyralis', 'Ignelater luminosus', 'Rhynchophorus ferrugineus', 'Tribolium madens', 'Dendroctonus ponderosae', 'Dendroctonus ponderosae', 'Zophobas morio', 'Tenebrio molitor', 'Cylas formicarius', 'Acanthoscelides obtectus']


#################################################
############# make lineage tree #################
#################################################

# code from this github repository: https://github.com/slaide/generate-taxonimic-lineage-tree/blob/main/main.py

if True:

    ncbi = NCBITaxa()

    # ANSI escape code for blue color
    BLUE = '\033[94m'
    # ANSI escape code to reset color
    RESET = '\033[0m'

    class TreeNode:
        def __init__(
            self,
            node_name:str,
            tax_id:str,
            init_as_root_with_species:list[str]=None
        ):
            self.name = node_name
            self.taxid = tax_id
            self.children = {}

            # Prepare the data structures
            lineages = {}
            names = {}

            species_names=init_as_root_with_species

            if species_names is not None:
                # Query NCBI and fill `lineages` and `names`
                for species in species_names:
                    taxid = ncbi.get_name_translator([species]).get(species, [None])[0]
                    if taxid:
                        lineage_ids = ncbi.get_lineage(taxid)
                        lineage_names = ncbi.get_taxid_translator(lineage_ids)
                        lineages[species] = lineage_ids
                        names.update(lineage_names)

                # Insert lineages into the tree
                for species, lineage in lineages.items():
                    self.insert_lineage(lineage, names)

        def add_child(self, child):
            if child.taxid not in self.children:
                self.children[child.taxid] = child
        
        def __repr__(self):
            return f"{self.name} ({self.taxid})"

        def insert_lineage(self, lineage, names):
            current_node = self
            for taxid in lineage[1:]:  # skip root as it's already created
                if taxid not in current_node.children:
                    new_node = TreeNode(names[taxid], taxid)
                    current_node.add_child(new_node)
                current_node = current_node.children[taxid]

        def print(self, depth=0):
            """
            print a simple representation of this tree

            the species names in the output are highlighted in color.
            """

            if len(self.children)==0:
                print(BLUE,end="")
            print("  " * depth + str(self))
            if len(self.children)==0:
                print(RESET,end="")

            for child in self.children.values():
                child.print(depth + 1)

        def print_compact_tree(self, depth=0, print_children=False):
            """
            Print a simplified version of the tree

            this shows only:
                - leafs
                - the most immediate common ancestor of any leaf or branch
                - i.e.: does NOT show intermediate branches

            the species names in the output are highlighted in color.
            """
            
            # Determine if the current node has leaf nodes as direct children
            has_species_as_direct_children = any(len(child.children) == 0 for child in self.children.values())
            
            # Construct the node label, highlighting species names in blue
            label = f"{self.name} ({self.taxid})"
            if len(self.children) == 0:  # This is a species node
                label = BLUE + label + RESET
            
            if print_children or has_species_as_direct_children:
                print("  " * depth + label)
            
            for child in self.children.values():
                # If the current node has species as direct children, print all children (species).
                # Otherwise, pass the flag down to print children only if they have species as direct children.
                child.print_compact_tree(depth + 1, print_children=has_species_as_direct_children or print_children)


        def build_ete_tree(self, parent=None)->Tree:
            """
            Recursively builds an ETE tree, attaching taxonomic information to each node.
            """
            # If this is the root node
            if parent is None:
                ete_node = Tree(name=self.name)
            else:
                ete_node = parent.add_child(name=self.name)

            # Attach a hypothetical taxonomic identifier to the node as an example of taxonomic information.
            # In a real scenario, this could be a NCBI Taxonomy ID or other relevant identifier.
            ete_node.add_features(taxid="TaxID_" + self.name.replace(" ", "_"))

            for child in self.children.values():
                child.build_ete_tree(parent=ete_node)

            return ete_node

    def insert_lineage(root, lineage, names):
        current_node = root
        for taxid in lineage[1:]:  # skip root as it's already created
            if taxid not in current_node.children:
                new_node = TreeNode(names[taxid], taxid)
                current_node.add_child(new_node)
            current_node = current_node.children[taxid]

    def print_tree(node, depth=0):
        print("  " * depth + str(node))
        for child in node.children.values():
            print_tree(child, depth + 1)


    # added the callosobruchuses and D. melanogaster
    # species_names = ['Megabruchidius dorsalis', 'Megabruchidius tonkineus', 'Bruchidius siliquastri', 'Dendroctonus ponderosae', 'Anoplophora glabripennis', 'Leptinotarsa decemlineata', 'Onthophagus taurus', 'Agrilus planipennis', 'Nicrophorus vespilloides', 'Sitophilus oryzae', 'Diabrotica virgifera', 'Asbolus verrucosus', 'Photinus pyralis', 'Ignelater luminosus', 'Rhynchophorus ferrugineus', 'Tribolium madens', 'Dendroctonus ponderosae', 'Dendroctonus ponderosae', 'Zophobas morio', 'Tenebrio molitor', 'Cylas formicarius', 'Acanthoscelides obtectus', 'Callosobruchus analis', 'Callosobruchus chinensis', 'Callosobruchus maculatus', 'Tribolium castaneum', 'Drosophila melanogaster']
    # species_names = ['Megabruchidius dorsalis', 'Megabruchidius tonkineus', 'Bruchidius siliquastri', 'Dendroctonus ponderosae', 'Diabrotica virgifera', 'Asbolus verrucosus', 'Photinus pyralis', 'Ignelater luminosus', 'Rhynchophorus ferrugineus', 'Dendroctonus ponderosae', 'Zophobas morio', 'Tenebrio molitor', 'Acanthoscelides obtectus', 'Callosobruchus analis', 'Callosobruchus chinensis', 'Callosobruchus maculatus', 'Tribolium castaneum', 'Drosophila melanogaster', 'Megabruchidius tonkineus', 'Megabruchidius dorsalis', 'Coccinella septempunctata']

if __name__ == "__main__":
    
    data = "/Users/miltr339/work/PhD_code/PhD_chapter1/data/"
    species_names = ["Acanthoscelides obtectus", "Asbolus verrucosus", "Bruchidius siliquastri", "Callosobruchus chinensis", "Callosobruchus maculatus", "Coccinella septempunctata", "Dendroctonus ponderosae", "Ignelater luminosus", "Photinus pyralis", "Rhynchophorus ferrugineus", "Tenebrio molitor", "Tribolium castaneum", "Zophobas morio"] 
    print(len(species_names))

    root = TreeNode("root", 1, species_names)

    # Print the tree
    print("basic (custom) tree representation:")
    root.print()
    # print compact version of the tree
    print("compact (custom) tree representation:")
    root.print_compact_tree()


    # construct ete3 tree from custom data structure
    ete_tree = root.build_ete_tree()

    # print ete3 tree
    print("basic newick tree:")
    print(ete_tree.write(format=1))

    # For visualization, print the tree showing taxid annotations
    # this is a quite useful representation, with more information than the default newick tree 
    #print("well-formatted ASCII ete3 tree:")
    #print(ete_tree.get_ascii(show_internal=True, attributes=["name", "taxid"]))

    # save the tree as newick tree to a file, and include the taxid feature in each leaf
    ete_tree.write(outfile=f"{data}lineage_tree_with_taxid.newick",features=["taxid", "path_score"])


#################################################
############# plot lineage tree #################
#################################################

# copy-pasted command line result from above

# (((((((((Acanthoscelides obtectus:1)Acanthoscelides:1,(Bruchidius siliquastri:1)Bruchidius:1,(Callosobruchus chinensis:1,Callosobruchus maculatus:1)Callosobruchus:1)Bruchini:1)Bruchinae:1)Chrysomelidae:1)Chrysomeloidea:1,((((Asbolus verrucosus:1)Asbolus:1)Pimeliinae:1,((Tenebrio molitor:1)Tenebrio:1,(Zophobas morio:1)Zophobas:1)Tenebrioninae:1,((Tribolium castaneum:1)Tribolium:1)Tenebrionidae incertae sedis:1)Tenebrionidae:1)Tenebrionoidea:1,(((((Coccinella septempunctata:1)Coccinella:1)Coccinellini:1)Coccinellinae:1)Coccinellidae:1)Coccinelloidea:1,((((Dendroctonus ponderosae:1)Dendroctonus:1)Scolytinae:1,((Rhynchophorus ferrugineus:1)Rhynchophorus:1)Dryophthorinae:1)Curculionidae:1)Curculionoidea:1)Cucujiformia:1,((((((Ignelater luminosus:1)Ignelater:1)Pyrophorini:1)Agrypninae:1)Elateridae:1,(((Photinus pyralis:1)Photinus:1)Lampyrinae:1)Lampyridae:1)Elateroidea:1)Elateriformia:1)Polyphaga:1)Coleoptera:1);
















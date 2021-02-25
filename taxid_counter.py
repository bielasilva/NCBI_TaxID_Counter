#!/usr/bin/env python3

from Bio import Entrez
import csv, os, sys

def main():
    """Function to  run the program."""
    print("\n## Beginning program ##")
    if sys.argv[1] != "-compare":
        print("Normal mode selected")
        for file in sys.argv[1:]:
            core_out = core(file)
            print(f"Done with {core_out[0]} \n")
    if sys.argv[1] == "-compare":
        print("Compare mode selected")
        compare()

def core(file):
    """Program core, runs the standard operations (getting the TaxIDs and counting them"""
    metag_result_file = str(os.path.abspath(file))
    file_name = metag_result_file.split("/")[-1]
    get_taxids_out = get_taxids(metag_result_file, file_name)
    get_tax_data_out = get_tax_data(get_taxids_out)
    sort_tax_data_out = sort_tax_data(get_tax_data_out, file_name)
    save_file(file_name, sort_tax_data_out)
    return file_name, sort_tax_data_out

def compare():
    print("Entering compare module")
    complete_data = {}
    for file in sys.argv[2:]:
        core_out = core(file)
        if "centrifuge" in core_out[0]:
            complete_data["centrifuge"] = core_out[1]
        if "krk" in core_out[0]:
            complete_data["kraken"] = core_out[1]
        if "metacache" in core_out[0]:
            complete_data["metacache"] = core_out[1]
        if "mpa-out" in core_out[0]:
            complete_data["metaphlan"] = core_out[1]
        print(f"Done with {core_out[0]} \n")
    # with open("complete_data.txt", "w+") as out:
    #     print(complete_data, file=out)
    shared_bacteria_taxid = {}
    shared_fungi_taxid = {}
    shared_viruses_taxid = {}
    shared_archaea_taxid = {}
    ranks = ["phylum", "class", "order", "family", "genus", "species"]
    # for rank in ranks:
    #     shared_bacteria_taxid[rank] = complete_data["centrifuge"]["bacteria"][rank].intersection(complete_data["kraken"]["bacteria"][rank],
    #                                                                                                 complete_data["metacache"]["bacteria"][rank],
    #                                                                                                 complete_data["metaphlan"]["bacteria"][rank])
    #     shared_archaea_taxid[rank] = complete_data["centrifuge"]["archaea"][rank].intersection(complete_data["kraken"]["archaea"][rank],
    #                                                                                                 complete_data["metacache"]["archaea"][rank],
    #                                                                                                 complete_data["metaphlan"]["archaea"][rank])
    #     shared_viruses_taxid[rank] = complete_data["centrifuge"]["viruses"][rank].intersection(complete_data["kraken"]["viruses"][rank],
    #                                                                                                 complete_data["metacache"]["viruses"][rank],
    #                                                                                                 complete_data["metaphlan"]["viruses"][rank])
    #     shared_fungi_taxid[rank] = complete_data["centrifuge"]["fungi"][rank].intersection(complete_data["kraken"]["fungi"][rank],
    #                                                                                                 complete_data["metacache"]["fungi"][rank],
    #                                                                                                 complete_data["metaphlan"]["fungi"][rank])
    for rank in ranks:
        shared_bacteria_taxid[rank] = complete_data["centrifuge"]["bacteria"][rank].intersection(complete_data["kraken"]["bacteria"][rank])
        shared_archaea_taxid[rank] = complete_data["centrifuge"]["archaea"][rank].intersection(complete_data["kraken"]["archaea"][rank])
        shared_viruses_taxid[rank] = complete_data["centrifuge"]["viruses"][rank].intersection(complete_data["kraken"]["viruses"][rank])
        shared_fungi_taxid[rank] = complete_data["centrifuge"]["fungi"][rank].intersection(complete_data["kraken"]["fungi"][rank])
    total = 0
    print("\nBacteria")
    for key, value in shared_bacteria_taxid.items():
        print(key, "->", len(value))
        total += len(value)
    print(f"Total = {total}")
    total = 0
    
    print("\nArchaea")
    for key, value in shared_archaea_taxid.items():
        print(key, "->", len(value))
        total += len(value)
    print(f"Total = {total}")
    
    total = 0
    print("\nFungi")
    for key, value in shared_fungi_taxid.items():
        print(key, "->", len(value))
        total += len(value)
    print(f"Total = {total}")
    
    total = 0
    print("\nViruses")
    for key, value in shared_viruses_taxid.items():
        print(key, "->", len(value))
        total += len(value)
    print(f"Total = {total}")

    # shared_bacteria_names = {}
    # shared_fungi_names = {}
    # shared_viruses_names = {}
    # shared_archaea_names = {}
    # for key, value in shared_bacteria_taxid.items():
    #     print(f"Getting data from Bacteria {key}")
    #     get_tax_data_out = get_tax_data(list(value))
    #     shared_bacteria_names[f'{key}'] = [item["ScientificName"] for item in get_tax_data_out]
    # for key, value in shared_archaea_taxid.items():
    #     print(f"Getting data from archaea {key}")
    #     get_tax_data_out = get_tax_data(list(value))
    #     shared_archaea_names[f'{key}'] = [item["ScientificName"] for item in get_tax_data_out]
    # for key, value in shared_fungi_taxid.items():
    #     print(f"Getting data from fungi {key}")
    #     get_tax_data_out = get_tax_data(list(value))
    #     shared_fungi_names[f'{key}'] = [item["ScientificName"] for item in get_tax_data_out]
    # for key, value in shared_viruses_taxid.items():
    #     print(f"Getting data from viruses {key}")
    #     get_tax_data_out = get_tax_data(list(value))
    #     shared_viruses_names[f'{key}'] = [item["ScientificName"] for item in get_tax_data_out]
        
    # print(shared_bacteria_names,"\n", shared_archaea_names,"\n", shared_fungi_names,"\n", shared_viruses_names)

# This function gets the TaxIDs from the report file
def get_taxids(in_file, origin):
    """Gets the output file from the classification program and returns a list with all the TaxIDs."""
    print(f"\nGetting TaxIDs from {origin}")
    taxid_list = set([])
    with open(in_file) as csv_read:
        mg_result = csv.reader(csv_read, delimiter='\t')
        metaphlan = False
        if any(st in in_file for st in ["centrifuge", "krk"]): #centrifuge, kraken2 and metaplhan3
            index = 2
        elif "metacache" in in_file: #metacache
            index = 3
        elif "mpa-out" in in_file:
            index = 1
            metaphlan = True
        else:
            print("File format not recognized")
            exit()
        for row in mg_result:
            if not row[0].startswith("#"):
                if metaphlan:
                    for item in row[index].split("|"):
                        taxid_list.add(item)
                elif row[index] != "0":
                    taxid_list.add(row[index])
        taxid_list = list(taxid_list)
    print(f"Got {len(taxid_list)} unique TaxIDs")
    return taxid_list

# This function gets the taxonomy data from the TaxIDs
def get_tax_data(list):
    """Takes a list with TaxIDs and returns a list with the data from NCBI Taxonomy."""
    print("Getting TaxIDs data")
    Entrez.email = 'gabrielamorimsilva@gmail.com'
    Entrez.api_key = 'f19b2580b4e240476bdef13bab28f8bb7808'
    taxid_data_list = []
    sub_taxid_list = [list[i:i + 10000] for i in range(0, len(list), 10000)]
    for sub_list in sub_taxid_list:
        sub_taxid_str = ','.join(sub_list)
        search = Entrez.efetch(id=sub_taxid_str, db="taxonomy", retmode="xml")
        taxid_data_list.extend(Entrez.read(search))
    print(f"Got data from {len(taxid_data_list)}")
    return taxid_data_list

# This function sorts and categorize the taxonomy data 
def sort_tax_data(data_list, origin):
    """Takes the data from the TaxIDs and sorts it to each taxonomy level and category (Bacteria, archaea, viruses, fungi and protozoa) and returns it."""
    print("Sorting Tax data")

    bacteria = {"superkingdom": set([]), "kingdom": set([]), "phylum": set([]), "class": set([]), "order": set([]), "family": set([]), "genus": set([]), "species": set([])}
    archaea = {"superkingdom": set([]), "kingdom": set([]), "phylum": set([]), "class": set([]), "order": set([]), "family": set([]), "genus": set([]), "species": set([])}
    viruses = {"superkingdom": set([]), "kingdom": set([]), "phylum": set([]), "class": set([]), "order": set([]), "family": set([]), "genus": set([]), "species": set([])}
    fungi = {"superkingdom": set([]), "kingdom": set([]), "phylum": set([]), "class": set([]), "order": set([]), "family": set([]), "genus": set([]), "species": set([])}
    protozoa = {"superkingdom": set([]), "kingdom": set([]), "phylum": set([]), "class": set([]), "order": set([]), "family": set([]), "genus": set([]), "species": set([])}
    rest = set([])
    ranks = ["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
    
    try:
        for entry in data_list:
            trying = entry
            if entry["TaxId"] in ("1", "131567", "10239", "2157"):
                if entry["TaxId"] in ("1", "131567"):
                    rest.add(entry["TaxId"])
                elif entry["TaxId"] == "10239":
                    viruses["superkingdom"].add(entry["TaxId"])
                elif entry["TaxId"] == "2157":
                    archaea["superkingdom"].add(entry["TaxId"])
            elif "2" in (entry["TaxId"], (next((item["TaxId"] for item in entry["LineageEx"] if item["TaxId"] == "2"), None))):
                for rank in ranks:
                    if (next((item["TaxId"] for item in entry["LineageEx"] if item["Rank"] == rank), None)) != None:
                        bacteria[rank].add((next((item["TaxId"] for item in entry["LineageEx"] if item["Rank"] == rank), None)))
            elif "2157" in (entry["TaxId"], (next((item["TaxId"] for item in entry["LineageEx"] if item["TaxId"] == "2157"), None))):
                for rank in ranks:
                    if (next((item["TaxId"] for item in entry["LineageEx"] if item["Rank"] == rank), None)) != None:
                        archaea[rank].add((next((item["TaxId"] for item in entry["LineageEx"] if item["Rank"] == rank), None)))
            elif "10239" in (entry["TaxId"], (next((item["TaxId"] for item in entry["LineageEx"] if item["TaxId"] == "10239"), None))):
                for rank in ranks:
                    if (next((item["TaxId"] for item in entry["LineageEx"] if item["Rank"] == rank), None)) != None:
                        viruses[rank].add((next((item["TaxId"] for item in entry["LineageEx"] if item["Rank"] == rank), None)))
            elif "4751" in (entry["TaxId"], (next((item["TaxId"] for item in entry["LineageEx"] if item["TaxId"] == "4751"), None))):
                for rank in ranks:
                    if (next((item["TaxId"] for item in entry["LineageEx"] if item["Rank"] == rank), None)) != None:
                        fungi[rank].add((next((item["TaxId"] for item in entry["LineageEx"] if item["Rank"] == rank), None)))
            elif ("5794" in (entry["TaxId"], (next((item["TaxId"] for item in entry["LineageEx"] if item["TaxId"] == "5794"), None)))
                    or "33682" in (entry["TaxId"], (next((item["TaxId"] for item in entry["LineageEx"] if item["TaxId"] == "33682"), None)))):
                for rank in ranks:
                    if (next((item["TaxId"] for item in entry["LineageEx"] if item["Rank"] == rank), None)) != None:
                        protozoa[rank].add((next((item["TaxId"] for item in entry["LineageEx"] if item["Rank"] == rank), None)))
            else:
                rest.add(entry["TaxId"])
    except:
        print(f"falty is {trying}")
    
    result_dict = {"bacteria": bacteria, "archaea": archaea, "viruses": viruses, "fungi": fungi, "protozoa": protozoa, "rest": rest}
    return result_dict

def save_file(file_name, result_dict):
    """Saves the data to a CSV file."""
    with open(f"{file_name}_sorted.csv", "w+") as csv_out:
        csv_writer = csv.writer(csv_out, delimiter=',')
        csv_writer.writerow(["", "Bacteria", "Archaea", "Viruses", "Fungi", "Protozoa"])
        ranks = ["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"]
        for rank in ranks:
            csv_writer.writerow([rank,
                                len(result_dict["bacteria"][rank]),
                                len(result_dict["archaea"][rank]),
                                len(result_dict["viruses"][rank]),
                                len(result_dict["fungi"][rank]),
                                len(result_dict["protozoa"][rank])])
        csv_writer.writerow(["Not classified", result_dict["rest"], "", "", "",""])
    # with open(f"{file_name}_sorted.out", "w+") as out:
    #     out.write(result_str)

if __name__ == '__main__':
    main()

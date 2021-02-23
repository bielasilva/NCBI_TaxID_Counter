#!/usr/bin/env python3

from Bio import Entrez
import csv, os, sys

def main():
    """Function to  run the program."""
    print("\n## Beginning program ##")
    for file in sys.argv[1:]:
        metag_result_file = str(os.path.abspath(file))
        file_name = metag_result_file.split("/")[-1]
        get_taxids_out = get_taxids(metag_result_file, file_name)
        get_tax_data_out = get_tax_data(get_taxids_out)
        sort_tax_data_out = sort_tax_data(get_tax_data_out, file_name)
        save_file(file_name, sort_tax_data_out[0], sort_tax_data_out[1])
        print(f"Done with {file_name} \n")


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
                    if (next((item["TaxId"] for item in entry["LineageEx"] if item["Rank"] == f"{rank}"), None)) != None:
                        bacteria[f"{rank}"].add((next((item["TaxId"] for item in entry["LineageEx"] if item["Rank"] == f"{rank}"), None)))
            elif "2157" in (entry["TaxId"], (next((item["TaxId"] for item in entry["LineageEx"] if item["TaxId"] == "2157"), None))):
                for rank in ranks:
                    if (next((item["TaxId"] for item in entry["LineageEx"] if item["Rank"] == f"{rank}"), None)) != None:
                        archaea[f"{rank}"].add((next((item["TaxId"] for item in entry["LineageEx"] if item["Rank"] == f"{rank}"), None)))
            elif "10239" in (entry["TaxId"], (next((item["TaxId"] for item in entry["LineageEx"] if item["TaxId"] == "10239"), None))):
                for rank in ranks:
                    if (next((item["TaxId"] for item in entry["LineageEx"] if item["Rank"] == f"{rank}"), None)) != None:
                        viruses[f"{rank}"].add((next((item["TaxId"] for item in entry["LineageEx"] if item["Rank"] == f"{rank}"), None)))
            elif "4751" in (entry["TaxId"], (next((item["TaxId"] for item in entry["LineageEx"] if item["TaxId"] == "4751"), None))):
                for rank in ranks:
                    if (next((item["TaxId"] for item in entry["LineageEx"] if item["Rank"] == f"{rank}"), None)) != None:
                        fungi[f"{rank}"].add((next((item["TaxId"] for item in entry["LineageEx"] if item["Rank"] == f"{rank}"), None)))
            elif ("5794" in (entry["TaxId"], (next((item["TaxId"] for item in entry["LineageEx"] if item["TaxId"] == "5794"), None)))
                    or "33682" in (entry["TaxId"], (next((item["TaxId"] for item in entry["LineageEx"] if item["TaxId"] == "33682"), None)))):
                for rank in ranks:
                    if (next((item["TaxId"] for item in entry["LineageEx"] if item["Rank"] == f"{rank}"), None)) != None:
                        protozoa[f"{rank}"].add((next((item["TaxId"] for item in entry["LineageEx"] if item["Rank"] == f"{rank}"), None)))
            else:
                rest.add(entry["TaxId"])
    except:
        print(f"falty is {trying}")
    
    result_dict = {"bacteria": bacteria, "archaea": archaea, "viruses": viruses, "fungi": fungi, "protozoa": protozoa, "rest": rest}
    result_str = f"""
    {origin}
    Bacteria has 
        Superkingdom: {len(bacteria["superkingdom"])}
        Kingdom: {len(bacteria["kingdom"])}
        Phylum: {len(bacteria["phylum"])}
        Class: {len(bacteria["class"])}
        Order: {len(bacteria["order"])}
        Family: {len(bacteria["family"])}
        Genus: {len(bacteria["genus"])}
        Species: {len(bacteria["species"])}

    Archaea has 
        Superkingdom: {len(archaea["superkingdom"])}
        Kingdom: {len(archaea["kingdom"])}
        Phylum: {len(archaea["phylum"])}
        Class: {len(archaea["class"])}
        Order: {len(archaea["order"])}
        Family: {len(archaea["family"])}
        Genus: {len(archaea["genus"])}
        Species: {len(archaea["species"])}

    Viruses has 
        Superkingdom: {len(viruses["superkingdom"])}
        Kingdom: {len(viruses["kingdom"])}
        Phylum: {len(viruses["phylum"])}
        Class: {len(viruses["class"])}
        Order: {len(viruses["order"])}
        Family: {len(viruses["family"])}
        Genus: {len(viruses["genus"])}
        Species: {len(viruses["species"])}

    Fungi has 
        Superkingdom: {len(fungi["superkingdom"])}
        Kingdom: {len(fungi["kingdom"])}
        Phylum: {len(fungi["phylum"])}
        Class: {len(fungi["class"])}
        Order: {len(fungi["order"])}
        Family: {len(fungi["family"])}
        Genus: {len(fungi["genus"])}
        Species: {len(fungi["species"])}

    Protozoa has 
        Superkingdom: {len(protozoa["superkingdom"])}
        Kingdom: {len(protozoa["kingdom"])}
        Phylum: {len(protozoa["phylum"])}
        Class: {len(protozoa["class"])}
        Order: {len(protozoa["order"])}
        Family: {len(protozoa["family"])}
        Genus: {len(protozoa["genus"])}
        Species: {len(protozoa["species"])}

    The {len(rest)} TaxIDs bellow were not classified:
        {rest}
    """
    
    return result_dict, result_str

def save_file(file_name, result_dict, result_str):
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

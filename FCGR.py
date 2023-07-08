import sys
import os                                                           
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from Bio import SeqIO


def fasta_reader(fasta):
	"""Return a generator with sequence description and sequence

	:param fasta: str filename
	"""
	flist = SeqIO.parse(fasta, "fasta")
	for i in flist:
		yield i.description, i.seq
nucleotides=["A","T","G","C","a","t","g","c"]

def get_kmer_frequencies(sequence, k):
    """
    Calculate k-mer frequencies for a given sequence.

    :param sequence: str, the DNA sequence to calculate k-mer frequencies for
    :param k: int, the length of k-mer to calculate frequencies for
    :return: dict, a dictionary of k-mers and their frequencies
    """
    freq = {}
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        if kmer not in freq:
            freq[kmer] = 0
        freq[kmer] += 1
    for key in list(freq):
        if not all(x in nucleotides for x in key):
            del freq[key]
    return freq

def generate_coordinates(k):
    """
    Generate coordinates for the Frequency Chaos Game Representation (FCGR) based on the length of the k-mer.
    Not really used to create the FCGR but may be useful to keep Cartesian coordinates.

    :param k: int, the length of k-mer
    :return: dict, a dictionary of coordinates for each nucleotide
    """
    coords = {"G": (0.5, 0.5), "C": (0.5, -0.5), "A": (-0.5, -0.5), "T": (-0.5, 0.5)}
    for i in range(k - 1):
        new_coords = {}
        for nuc, (x, y) in coords.items():
            new_coords[nuc + "A"] = ((x / 2) + 0.5, (y / 2) + 0.5)
            new_coords[nuc + "T"] = ((x / 2) + 0.5, (y / 2) - 0.5)
            new_coords[nuc + "G"] = ((x / 2) - 0.5, (y / 2) - 0.5)
            new_coords[nuc + "C"] = ((x / 2) - 0.5, (y / 2) + 0.5)
        coords.update(new_coords)
    fix_coords = {i:coords[i] for i in coords if len(i)>=k}
    coords.update(fix_coords)
    return fix_coords
 
def chaos_game_representation(frequencies, k):
    """Generate the Frequency Chaos Game Representation of a given DNA sequence.

    Args:
        frequencies (dict): Dictionary containing the Kmers and their frequencies
        k (int): k-mer length

    Returns:
        np.ndarray: A 2D numpy array with the FCGR.
    """
    fcgr = np.zeros((2**k, 2**k))
    maxx = maxy = 2**k
    posx = 0
    posy = 0
    max_count = max(frequencies.values())
    for key, value in frequencies.items():
        for char in key:
            if char == "T":
                posx += maxx // 2
            elif char == "C":
                posy += maxy // 2
            elif char == "G":
                posx += maxx // 2
                posy += maxy // 2
            maxx = maxx // 2
            maxy = maxy // 2
        fcgr[(posx),(posy)] = value / max_count
        maxx = maxy = 2**k
        posx = 0
        posy = 0
    return fcgr

def get_args():
    """Get args"""
    import argparse
    parser = argparse.ArgumentParser(
        description="Chaos Game Representation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
        )
    parser.add_argument("files", nargs="*")
    parser.add_argument("-k", "--kmers", action="store", default="3",help="Determine kmers, default is 3") 
    parser.add_argument("-d", "--dpi" , action="store_true", default="1", help="set resolution of FCGR image")  
    parser.add_argument("-g", "--generate_plot", action="store_true", help="Generate plot if this flag is set")
    args = parser.parse_args()
    if not args.files:
        raise TypeError("Not Enough Argument; expected at least 1")
    for i in args.files:
        if not os.path.isfile(i):
            raise TypeError(
                "File %s does not exist or is not the right type" % i
                )
    return args


if __name__ == "__main__":
    fig_id = 1
    args = get_args()
    for i in args.files:
        original_file_name = os.path.basename(i)
        fasta_seq = fasta_reader(i)
        concatenated_seq = ""
        for name, seq in fasta_seq:
            print("Procesando secuencia:",name)
            name = name.replace(":","(colon)") #Necessary if they use a name that denotes coordinates.
            concatenated_seq += str(seq)
        if args.kmers:
            k = int(args.kmers)
        if args.dpi:
            user_res = int(args.dpi)
        freq = get_kmer_frequencies(concatenated_seq, k)
        print(freq)
        # generate coordinates if you want
        #coords = generate_coordinates(k)
        #print("Generated coords", coords)

        chaos_k4 = chaos_game_representation(freq, k)
        print(chaos_k4)
        print(chaos_k4.shape)
        flattened_arr= chaos_k4.flatten()
        print(flattened_arr)
        #
        np.savetxt(original_file_name+"_FCGR" + str(k) + "-mers_preprocessed.txt",  flattened_arr, newline=",", fmt="%.16f")
        if args.generate_plot:
            plt.title("Frequency chaos game representation for " + str(k) + "-mers")
            # Add letters outside each corner of the plot
            plt.text(-0.1, -0.1, "A", ha="right", va="top", transform=plt.gca().transAxes, fontsize=14)
            plt.text(1.1, -0.1, "C", ha="left", va="top", transform=plt.gca().transAxes, fontsize=14)
            plt.text(-0.1, 1.1, "T", ha="right", va="bottom", transform=plt.gca().transAxes, fontsize=14)
            plt.text(1.1, 1.1, "G", ha="left", va="bottom", transform=plt.gca().transAxes, fontsize=14)
            plt.imshow(chaos_k4, interpolation="nearest", cmap=cm.gray_r, origin="lower")
            plt.savefig(original_file_name+str(-k)+"-mers.png",bbox_inches="tight", pad_inches=0)
            fig = plt.figure(figsize=(chaos_k4.shape[1], chaos_k4.shape[0]), dpi=user_res)
            ax = plt.Axes(fig, [0., 0., 1., 1.], )
            ax.set_axis_off()
            fig.add_axes(ax)
            ax.imshow(chaos_k4, cmap=cm.gray_r, origin="lower")

            # Save the figure as a grayscale image without any extra elements
            plt.savefig(original_file_name+str(k)+"-mers_grayscale_image.png", dpi=user_res, bbox_inches="tight", pad_inches=0)
            plt.close(fig)


        
        


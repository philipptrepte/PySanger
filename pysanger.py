import os 
import sys 
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio.Seq import Seq
from snapgene_reader import snapgene_file_to_dict
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd 
import numpy as np
#import logomaker

matplotlib.rcParams['font.family']       = 'sans-serif'
matplotlib.rcParams['font.sans-serif']   = ["Arial","DejaVu Sans","Lucida Grande","Verdana"]
matplotlib.rcParams['figure.figsize']    = [3,3]
matplotlib.rcParams['font.size']         = 10
matplotlib.rcParams["axes.labelcolor"]   = "#000000"
matplotlib.rcParams["axes.linewidth"]    = 1.0 
matplotlib.rcParams["xtick.major.width"] = 1.0
matplotlib.rcParams["ytick.major.width"] = 1.0

_atgc_dict = {0:"A", 1:"T", 2:"G", 3:"C"}

def abi_to_dict(filename):
    record   = SeqIO.read(filename,'abi')
    abi_data = {"conf":[],
                "channel":{"A":[],
                           "T":[],
                           "G":[],
                           "C":[],
                          },
                "_channel":{"A":[],
                            "T":[],
                            "G":[],
                            "C":[],
                          }
                }
    for i, (pos, conf) in enumerate(zip(record.annotations['abif_raw']["PLOC1"], record.annotations['abif_raw']["PCON1"])):
        if pos > 4 and pos < len(record.annotations['abif_raw']["DATA9"])-5: 
            abi_data["conf"].append(conf)
            abi_data["channel"]["G"].append(record.annotations['abif_raw']["DATA9"][pos])
            abi_data["channel"]["A"].append(record.annotations['abif_raw']["DATA10"][pos])
            abi_data["channel"]["T"].append(record.annotations['abif_raw']["DATA11"][pos])
            abi_data["channel"]["C"].append(record.annotations['abif_raw']["DATA12"][pos])

            abi_data["_channel"]["G"].append(record.annotations['abif_raw']["DATA9"][pos-5])
            abi_data["_channel"]["G"].append(record.annotations['abif_raw']["DATA9"][pos-3])
            abi_data["_channel"]["G"].append(record.annotations['abif_raw']["DATA9"][pos])
            abi_data["_channel"]["G"].append(record.annotations['abif_raw']["DATA9"][pos+3])
            abi_data["_channel"]["G"].append(record.annotations['abif_raw']["DATA9"][pos+5])
            
            abi_data["_channel"]["A"].append(record.annotations['abif_raw']["DATA10"][pos-5])
            abi_data["_channel"]["A"].append(record.annotations['abif_raw']["DATA10"][pos-3])
            abi_data["_channel"]["A"].append(record.annotations['abif_raw']["DATA10"][pos])
            abi_data["_channel"]["A"].append(record.annotations['abif_raw']["DATA10"][pos+3])
            abi_data["_channel"]["A"].append(record.annotations['abif_raw']["DATA10"][pos+5])
            
            abi_data["_channel"]["T"].append(record.annotations['abif_raw']["DATA11"][pos-5])
            abi_data["_channel"]["T"].append(record.annotations['abif_raw']["DATA11"][pos-3])
            abi_data["_channel"]["T"].append(record.annotations['abif_raw']["DATA11"][pos])
            abi_data["_channel"]["T"].append(record.annotations['abif_raw']["DATA11"][pos+3])
            abi_data["_channel"]["T"].append(record.annotations['abif_raw']["DATA11"][pos+5])

            abi_data["_channel"]["C"].append(record.annotations['abif_raw']["DATA12"][pos-5])
            abi_data["_channel"]["C"].append(record.annotations['abif_raw']["DATA12"][pos-3])
            abi_data["_channel"]["C"].append(record.annotations['abif_raw']["DATA12"][pos])
            abi_data["_channel"]["C"].append(record.annotations['abif_raw']["DATA12"][pos+3])
            abi_data["_channel"]["C"].append(record.annotations['abif_raw']["DATA12"][pos+5])

    return abi_data, filename

def generate_consensusseq(abidata):
    consensus_seq = "" 
    
    for values in zip(abidata["channel"]["A"], abidata["channel"]["T"], abidata["channel"]["G"], abidata["channel"]["C"]):
        consensus_seq += _atgc_dict[values.index(max(values))]
     
    return (consensus_seq, consensus_seq.translate(str.maketrans("ATGC","TACG"))[::-1]) 

def generate_pwm(abidata):
    _abidata, seq_filename = abidata
    pwm = {"A":[], "T":[], "G":[], "C":[]} 
    for values in zip(_abidata["channel"]["A"], _abidata["channel"]["T"], _abidata["channel"]["G"], _abidata["channel"]["C"]):
        v = 100000 / (sum(values)+1) 
        new_values = (v*values[0], v*values[1], v*values[2], v*values[3])
        new_values = list(map(int, new_values))
        
        while sum(new_values) < 100000:
            for i in range(len(new_values)):
                new_values[i] += 1
                if sum(new_values) == 100000:
                    break 
        
        pwm["A"].append(new_values[0])
        pwm["T"].append(new_values[1])
        pwm["G"].append(new_values[2])
        pwm["C"].append(new_values[3])
    
    pwm=pd.DataFrame(pwm)
    return pwm 

def _colorbar(ax, ref, matches=None, char=True, fontsize=2):
    bars = ax.bar(list(range(len(ref))), [0.9] * (len(ref)), width=1.0, edgecolor="#BBBBBB", linewidth=0.5, align="edge",bottom=0.05)
    ax.set_xlim(0,len(ref))
    ax.set_ylim(0,1.00)
    p = 0
    if matches is None:
        for bar, c in zip(bars,ref):
            bar.set_facecolor("w")
            if char == True:
                ax.text(p+0.5,0.45,c,va="center",ha="center",fontsize=fontsize,zorder=100)     
            p += 1
    
    else:
        for m, bar, c in zip(matches, bars, ref):
            if m == -1:
                bar.set_alpha(0.5)
                bar.set_facecolor("r")
                bar.set_edgecolor("#BBBBBB")
                bar.set_linewidth(0.5)
            else:
                bar.set_facecolor("w")
                bar.set_edgecolor("#BBBBBB")
                bar.set_linewidth(0.5)
            
            if char == True:
                ax.text(p+0.5,0.45,c,va="center",ha="center",fontsize=fontsize,zorder=100)     
            p += 1
    ax.set_xticks([])
    ax.set_yticks([])
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)
    #ax.patch.set_alpha(0.0)
    return ax

def alignment(abidata=None, template=None, strand=1):  
    if abidata is not None:
        _abidata, seq_filename = abi_to_dict(abidata)
    else:
        raise ValueError("Invalid ABI file. Please provide a valid ABI file.")
    
    avalues = _abidata["_channel"]["A"]
    tvalues = _abidata["_channel"]["T"]
    gvalues = _abidata["_channel"]["G"]
    cvalues = _abidata["_channel"]["C"]
    consensus_seq_set = generate_consensusseq(_abidata)
    
    if strand == 1: 
        subject = consensus_seq_set[0]
    
    if strand == -1:
        subject = consensus_seq_set[1] 
        avalues, tvalues, gvalues, cvalues = tuple(map(lambda x:list(reversed(x)), tvalues, avalues, cvalues, gvalues)) 
    
    if isinstance(template, str) and template.endswith((".gb", ".dna", ".xdna")):
        print(f"Reading vector map {template} as template.")
        if template.endswith(".gb"):
            """Reads a .gb file and transforms all sequences to upper case."""
            with open(template, "r") as file:
                sequence = SeqIO.read(file, "genbank")
                sequence.seq = sequence.seq.upper()
                name = str(template)
                #return sequence.seq, name
                
        elif template.endswith(".dna"):
            """Reads a .dna file and transforms all sequences to upper case."""
            data = snapgene_file_to_dict(filepath=template)
            sequence = Seq(data["seq"])
            sequence.seq = sequence.upper()
            name = str(template)
            #return sequence.seq, name
            
        elif template.endswith(".xdna"):
            """Reads a .xdna file and transforms all sequences to upper case."""
            with open(template, "rb") as file:
                sequence = next(XdnaIterator(file))
                sequence.seq = sequence.seq.upper()
                name = str(template)
                #return sequence.seq, name    
                
        else:
            raise ValueError("Invalid file format. Please provide a .gb, .dna, or .xdna file.")
    
    elif type(template) == str:
        print("Reading user-provided sequence as template.")
        sequence = template.upper()
        name = "User-provided sequence"

    else:
        raise ValueError("Invalid template. Please provide a DNA string ('ATGC') or a file path.") 

    if sequence is not None:
        if name == "User-provided sequence":
            print(f"Aligning sequencing file `{seq_filename}` ({len(subject)} bp) to `{name}` ({len(sequence)} bp).")
            alignments  = pairwise2.align.globalms(sequence, subject, 2, 0, -10, -1, penalize_end_gaps=False)
        
        elif name != "User-provided sequence":
            print(f"Aligning sequencing file `{seq_filename}` ({len(subject)} bp) to provided reference sequence `{name}` ({len(sequence)} bp).")
            alignments = pairwise2.align.globalms(sequence, subject, 2, 0, -10, -1, penalize_end_gaps=False)

        atemplate   = alignments[0][0].upper()
        asubject    = alignments[0][1].upper()

        ts = 0 
        ss = 0 
        for t, s in zip(atemplate, asubject):
            if t != "-" and s == "-": 
                ts += 1
            
            elif t == "-" and s != "-": 
                ss += 1  
            
            elif t == "-" and s == "-":
                ts += 1
                ss += 1

            else:
                break 

        te = 0 
        se = 0 
        for t, s in zip(atemplate[::-1], asubject[::-1]):
            if t != "-" and s == "-": 
                te += 1
            elif t == "-" and s != "-": 
                se += 1  
            elif t == "-" and s == "-":
                te += 1
                se += 1
            else:
                break 
       
        pos = 0 
        tpos = 0
        matches = []
        new_avalues = []
        new_tvalues = [] 
        new_gvalues = []
        new_cvalues = [] 
        
        for t, s in zip(atemplate, asubject):
            if s == "-":
                new_avalues.extend([0,0,0,0,0]) 
                new_tvalues.extend([0,0,0,0,0]) 
                new_gvalues.extend([0,0,0,0,0]) 
                new_cvalues.extend([0,0,0,0,0])     
            else:
                new_avalues.append(avalues[5*pos]) 
                new_avalues.append(avalues[5*pos+1]) 
                new_avalues.append(avalues[5*pos+2])
                new_avalues.append(avalues[5*pos+3]) 
                new_avalues.append(avalues[5*pos+4])

                new_tvalues.append(tvalues[5*pos])
                new_tvalues.append(tvalues[5*pos+1]) 
                new_tvalues.append(tvalues[5*pos+2])
                new_tvalues.append(tvalues[5*pos+3]) 
                new_tvalues.append(tvalues[5*pos+4])
                
                new_gvalues.append(gvalues[5*pos])
                new_gvalues.append(gvalues[5*pos+1]) 
                new_gvalues.append(gvalues[5*pos+2])
                new_gvalues.append(gvalues[5*pos+3]) 
                new_gvalues.append(gvalues[5*pos+4])

                new_cvalues.append(cvalues[5*pos])
                new_cvalues.append(cvalues[5*pos+1])
                new_cvalues.append(cvalues[5*pos+2])
                new_cvalues.append(cvalues[5*pos+3]) 
                new_cvalues.append(cvalues[5*pos+4])

                pos += 1
            tpos += 1

            if t == s: 
                matches.append(1) 
            elif len(subject) > len(sequence) and s != t and ss < tpos < ss+len(sequence):
                matches.append(-1)
            elif len(sequence) > len(subject) and s != t and ts < tpos < ts+len(subject):
                matches.append(-1)
            else:
                matches.append(0) 
            
            match_count = matches.count(1)
            mismatch_count = matches.count(-1)
            not_aligned_count = matches.count(0)
        print(f"{match_count} matches found between base pairs {ts} and {ts+len(subject)} in '{name}'. {mismatch_count} mismatches found. {not_aligned_count} base pairs not aligned.")
        
        avalues = new_avalues 
        tvalues = new_tvalues
        gvalues = new_gvalues
        cvalues = new_cvalues
    return(_abidata, asubject, subject, atemplate, sequence, matches, ss, se, ts, te, avalues, tvalues, gvalues, cvalues, name, seq_filename)

def adjust_quality_scores(tasubject, quality_scores):
    # Identify the positions of the insertions in the alignment
    insertion_positions = [i for i, base in enumerate(tasubject) if base == '-']
    
    # Adjust the quality scores by inserting placeholder values at the insertion positions
    adjusted_quality_scores = list(quality_scores)
    for pos in insertion_positions:
        adjusted_quality_scores.insert(pos, -1)  # Insert None or a specific low-quality score
    
    return adjusted_quality_scores

def visualize(alignment, region="aligned", fontsize=2):
    _abidata, asubject, subject, atemplate, sequence, matches, ss, se, ts, te, avalues, tvalues, gvalues, cvalues, name, seq_filename = alignment
    plt.close("all")

    if region == "all":
        tasubject  = asubject
        tatemplate = atemplate
        fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(25, 3), gridspec_kw={'height_ratios': [1, 3, 1]})

    elif region == "aligned":
        if(min(len(sequence), len(subject)) < 100):
            fig_width = 5
        else:
            fig_width = 20
        if ts == 0 and te == 0:
            tasubject  = asubject[ss:len(asubject) - se]
            tatemplate = atemplate[ss:len(atemplate) - se]
            matches   = matches[ss:len(asubject) - se]
        elif ts != 0 and te != 0:
            tasubject  = asubject[ts:len(asubject) - te]
            tatemplate = atemplate[ts:len(atemplate) - te]
            matches   = matches[ts:len(asubject) - te]
        fig, axs = plt.subplots(nrows=3, ncols=1, figsize=(fig_width, 3), gridspec_kw={'height_ratios': [1, 3, 1]})

    else:
        raise ValueError("Invalid region. Please provide 'all' or 'aligned'.")
    
    ax = axs[0]
    ax2 = axs[1]
    ax3 = axs[2]
    
    # Determine tick spacing based on the length of the sequences
    if region == "all":
        tlength = len(sequence)
    elif region == "aligned":
        tlength = len(tatemplate)
    if tlength < 100:
        ttick_space = 10
    elif tlength < 200:
        ttick_space = 20
    elif tlength < 500:
        ttick_space = 50
    elif tlength < 1000:
        ttick_space = 50
    elif tlength < 2500:
        ttick_space = 100    
    else:
        ttick_space = 250
    if ttick_space > len(sequence):
        if len(sequence) < 100:
            ttick_space = 10
        elif len(sequence) < 200:
            ttick_space = 20
        elif len(sequence) < 500:
            ttick_space = 50
        elif len(sequence) < 1000:
            ttick_space = 50
        elif len(sequence) < 2500:
            ttick_space = 100    

    if tatemplate is not None:
        ax = _colorbar(ax, tatemplate, matches=matches, char=True, fontsize=fontsize)
        tnum = 0
        if region == "all":
            if ts != 0:
                ss = ts
                ts = 0
            elif ts == 0:
                ts = ss
                ss = 0
            tticks      = [] 
            tticklabels = [] 
            for letter in tatemplate:
                if (tnum) % ttick_space == 0: 
                    tticks.append(ts + tnum) 
                    tticklabels.append(str(tnum))
                    tnum += 1
                if letter == "-":
                    pass 
                else:
                    tnum += 1 
        elif region == "aligned":
            tticks      = [0.5]
            tticklabels = [ts]
            for letter in tatemplate:
                if (tnum+1) % ttick_space == 0: 
                    tticks.append(tnum+0.5) 
                    tticklabels.append(str(ts+tnum+1))
                    tnum += 1
                if letter == "-":
                    pass 
                else:
                    tnum += 1 
        ax.xaxis.tick_top() 
        ax.set_xticks(tticks) 
        ax.set_xticklabels(tticklabels) 
        
        ax.set_ylabel(name, rotation=0, va="center", ha="right") 
    else:
        raise ValueError("Invalid template. No template sequence provided. Please check alignment.")
    
    if region == "all":
        quality = adjust_quality_scores(tasubject, quality_scores = _abidata["conf"])
    elif region == "aligned":
        quality = adjust_quality_scores(tasubject, quality_scores = _abidata["conf"])
    
    if region == "all":
        sns.barplot(x=list(range(len(tasubject))), y=quality, color="#F7F7FF", edgecolor="#BBBBBB", linewidth=0.5, ax=ax2, width = 1.0)
        ax2.set_xlim(-0.5, len(tatemplate) - 0.5) 
    elif region == "aligned":
        if ts == 0 and te == 0:
            sns.barplot(x=list(range(len(tasubject))), y=quality[ss:-se], color="#F7F7FF", edgecolor="#BBBBBB", linewidth=0.5, ax=ax2, width = 1.0)
        elif ts != 0 and te != 0:
            sns.barplot(x=list(range(len(tasubject))), y=quality, color="#F7F7FF", edgecolor="#BBBBBB", linewidth=0.5, ax=ax2, width = 1.0)
        ax2.set_xlim(-0.5, len(tasubject) - 0.5) 
    
    ax2X2 = ax2.twinx()

    if region == "all":
        min_length = min(len(tvalues), len(avalues), len(gvalues), len(cvalues))
        positions = np.linspace(-0.5, (min_length / 5) - 0.5, min_length).tolist()
        sns.lineplot(x=positions, y=tvalues, color="#FC58FE", lw=1, ax=ax2X2)
        sns.lineplot(x=positions, y=avalues, color="#33CC33", lw=1, ax=ax2X2)
        sns.lineplot(x=positions, y=gvalues, color="#303030", lw=1, ax=ax2X2)
        sns.lineplot(x=positions, y=cvalues, color="#395CC5", lw=1, ax=ax2X2)
        combined_values = tvalues + avalues + gvalues + cvalues

    elif region == "aligned":
        if ts == 0 and te == 0:
            _te = -5 * se
            ts = ss
        elif ts !=0 and te != 0:
            _te = -5 * te
        # Ensure the lengths of the arrays are consistent
        min_length = min(len(tvalues[5 * ts:_te]), len(avalues[5 * ts:_te]), len(gvalues[5 * ts:_te]), len(cvalues[5 * ts:_te]))
        positions = np.linspace(-0.5, (min_length / 5) - 0.5, min_length).tolist()
        sns.lineplot(x=positions, y=tvalues[5 * ts:5 * ts + min_length], color="#FC58FE", lw=1, ax=ax2X2)
        sns.lineplot(x=positions, y=avalues[5 * ts:5 * ts + min_length], color="#33CC33", lw=1, ax=ax2X2)
        sns.lineplot(x=positions, y=gvalues[5 * ts:5 * ts + min_length], color="#303030", lw=1, ax=ax2X2)
        sns.lineplot(x=positions, y=cvalues[5 * ts:5 * ts + min_length], color="#395CC5", lw=1, ax=ax2X2)

        combined_values = tvalues[5 * ss:5 * ss + min_length] + avalues[5 * ss:5 * ss + min_length] + gvalues[5 * ss:5 * ss + min_length] + cvalues[5 * ss:5 * ss + min_length]

    ax2.set_xticklabels("")
    ax2.set_xticks([])
    ax2.set_ylim(0, 1.01 * max(quality))
    ax2.set_ylabel("Quality")
    ax2X2.set_ylim(min(combined_values), 1.01 * max(combined_values))
    ax2X2.set_ylabel("Peak height")

    ax2.spines["top"].set_visible(False)
    ax2.spines["bottom"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2X2.spines["top"].set_visible(False)
    ax2X2.spines["bottom"].set_visible(False)
    ax2X2.spines["left"].set_visible(False)

    ax3 = _colorbar(ax3, tasubject, matches=matches, char=True, fontsize=fontsize)

    if region == "all":
        slength = len(subject)
    elif region == "aligned":
        slength = len(tasubject)
    if slength < 100:
        stick_space = 10
    elif slength < 200:
        stick_space = 20
    elif slength < 500:
        stick_space = 50
    elif slength < 1000:
        stick_space = 100
    elif slength < 2500:
        stick_space = 250    
    else:
        stick_space = 500
    snum = 0 
    if region == "all":
        sticks      = [ss] 
        sticklabels = [ss] 
        for letter in tasubject:
            if (ss+snum+1) % stick_space == 0: 
                sticks.append(ss+snum+1) 
                sticklabels.append(str(ss+snum+1))
                snum += 1
            if letter == "-":
                pass 
            else:
                snum += 1 
    elif region == 'aligned':
        if ts != 0 and te != 0:
            ts = ss
        sticks      = [0.5] 
        sticklabels = [ss]
        for letter in tasubject:
            if (ts + snum+1) % stick_space == 0: 
                sticks.append(snum+0.5) 
                sticklabels.append(str(ss+snum+1))
            if letter == "-":
                pass 
            else:
                snum += 1
    

    ax3.set_xticks(sticks) 
    ax3.set_xticklabels(sticklabels) 
    ax3.set_ylabel(seq_filename, rotation=0, va="top", ha="right") 

    # Adjust the spacing between subplots
    plt.subplots_adjust(hspace=0, wspace=0.3)  # Adjust hspace and wspace as needed

    #plt.tight_layout() # Adjust layout to ensure subplots fit into the figure area
    return fig

if __name__ == "__main__":

    align2        = alignment(abidata="BE MAFB5.ab1", template="AGCCGGCTGGCTGCAGGCGT")
    fig2          = visualize(align2, region = "all", fontsize = 5)
    fig2.show()
    fig2          = visualize(align2, region = "aligned", fontsize = 5)
    fig2.show()

    align        = alignment(abidata="seq_results/QPSQ0664-CMV-for.ab1", template="templates/QPPL0052_pcDNA3.1_mCitrine-C1-GW.dna")
    fig          = visualize(align, region="all")
    fig.show()
    fig          = visualize(align, region="aligned")
    fig.show()
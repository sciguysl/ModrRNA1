#!/usr/bin/env python
# coding: utf-8


import pandas as pd
import streamlit as st
import hashgen
import sys
sys.path.insert(0, './ncbi_taxonomy/')

st.set_page_config(
    page_title="ModrRNA1 - Bacterial rRNA Mod Finder",
    page_icon=":dna:",
#    layout="wide",

)

######################
# Page Title
######################
css='''
<style>
    section.main > div {max-width:75rem}
</style>
'''
st.markdown(css, unsafe_allow_html=True)



hide_streamlit_style = """
            <style>
            #MainMenu {visibility: hidden;}
            footer {visibility: hidden;}
            footer:after {
	content:'Made with love by Uvindu (BioTechie)'; 
	visibility: visible;
	display: block;
	position: relative;
	padding: 5px;
	top: 2px;
}
            </style>
            """
st.markdown(hide_streamlit_style, unsafe_allow_html=True)

st.markdown("""
        <style>
               .block-container {
                    padding-top: 1rem;
                    padding-bottom: 0rem;
                    padding-left: 2rem;
                    padding-right: 2rem;
                }
        </style>
        """, unsafe_allow_html=True)

#image = Image.open('dna-logo.jpg')
#

st.image("https://i.ibb.co/WtKyPnk/Untitled-design-1.png", use_column_width=True)
button_html = f'<a style="position: absolute; top: 0; right: 0;" href="mailto: uviwick@gmail.com" target="_blank">Report Errors</a>'
st.markdown(button_html, unsafe_allow_html=True)

st.write("""
# Bacterial rRNA Modification Sites Finder

This app detects the possible rRNA modification sites of any bacteria

***
""")


######################
# Input Text Box
######################

#st.sidebar.header('Enter DNA sequence')
st.header('Enter rRNA sequence')

sequence_input = ""
#sequence_input = " "

with st.form(key='my_form'):
    sequence = st.text_area("Sequence input", height=250, placeholder="Enter the rRNA sequence here")
    sequence = sequence.splitlines()
    sequence = ''.join(sequence) # Concatenates list to string
    submit_button = st.form_submit_button(label='Submit')
    
#sequence = st.sidebar.text_area("Sequence input", sequence_input, height=250)
#st.code(sequence)


st.write("""
***
""")

df_mods = pd.read_csv('bact_AllMods-processed.csv', names=["Chromosome", "Start", "End", "modID", "Score", "Strand", "Mod", "supportNum", "supportList", "supportListSub", "pubList", "cellList", "seqTypeList", "geneID", "transcriptID", "geneName", "geneType", "Region", "Seq", "Enzyme"])

df_mods = df_mods.iloc[1:]
df_mods['ID'] = None

for i in df_mods.index:
    input_string = df_mods['End'][i] + df_mods['modID'][i]
    unique_id = hashgen.generate_unique_id(input_string, length=3)
    df_mods['ID'][i] = unique_id


df_mods_fow = df_mods
df_mods_rev = df_mods.copy()


from Bio.Seq import Seq

for i in df_mods_rev.index:
  df_mods_rev['Seq'][i] = str(Seq(df_mods_rev['Seq'][i]).reverse_complement())


# In[12]


df_mods = pd.concat([df_mods_fow, df_mods_rev])
df_mods.reset_index(drop=True, inplace=True)


#

# In[18]:

refseq = sequence.upper()

if submit_button:
#    st.cache_data.clear()
    if sequence=="":
        st.error("Please enter the sequence")
        sys.exit()
    
#    sys.exit()

from Bio import pairwise2
from Bio.pairwise2 import format_alignment

def format_alignment2(align1, align2, score, begin, end, full_sequences=False):
    align_begin = begin
    align_end = end
    mod_position = 0
    start1 = start2 = ""
    start_m = begin  # Begin of match line (how many spaces to include)
    # For local alignments:
    if not full_sequences and (begin != 0 or end != len(align1)):
        # Calculate the actual start positions in the un-aligned sequences
        # This will only work if the gap symbol is '-' or ['-']!
        start1 = str(len(align1[:begin]) - align1[:begin].count("-") + 1) + " "
        start2 = str(len(align2[:begin]) - align2[:begin].count("-") + 1) + " "
        start_m = max(len(start1), len(start2))
    elif full_sequences:
        start_m = 0
        begin = 0
        end = len(align1)

    if isinstance(align1, list):
        # List elements will be separated by spaces, since they can be
        # of different lengths
        align1 = [a + " " for a in align1]
        align2 = [a + " " for a in align2]

    s1_line = ["{:>{width}}".format(start1, width=start_m)]  # seq1 line
    m_line = [" " * start_m]  # match line
    s2_line = ["{:>{width}}".format(start2, width=start_m)]  # seq2 line
    s3_line = ["{:>{width}}".format(start2, width=start_m)]
    linecounts=0
    for n, (a, b) in enumerate(zip(align1[begin:end], align2[begin:end])):
        # Since list elements can be of different length, we center them,
        # using the maximum length of the two compared elements as width

        if (b.strip() == "-"):
          linecounts = linecounts + 1
        startnum = len(align2[:begin]) - align2[:begin].count("-")
        m_len = max(len(a), len(b))
        s1_line.append("{:^{width}}".format(a, width=m_len))
        s2_line.append("{:^{width}}".format(b, width=m_len))
        s3_line.append("{:^{width}}".format(b, width=m_len))
        if full_sequences and (n < align_begin or n >= align_end):
            m_line.append("{:^{width}}".format(" ", width=m_len))  # space
            continue
        if (a == b and int(start2)+ n-linecounts == 21):
            mod_position = int(start1) + n-linecounts
            m_line.append("{:^{width}}".format("#", width=m_len))  # match
        elif (a == b):
            m_line.append("{:^{width}}".format("|", width=m_len))  # match
        elif a.strip() == "-" or b.strip() == "-":
            m_line.append("{:^{width}}".format(" ", width=m_len))  # gap
        else:
            m_line.append("{:^{width}}".format(".", width=m_len))  # mismatch
    s2_line.append(f"\n  Score={score:g}\n")
#    return "\n".join(["".join(s1_line), "".join(m_line), "".join(s2_line)])
    return mod_position


seqs = df_mods['Seq']

aligner = pairwise2.align.localms
high_conf_df = pd.DataFrame(columns=['Align', 'Desc', 'Score', 'Start', 'Enzyme', 'Mod', 'Mod2', 'ID'])
low_conf_df = pd.DataFrame(columns=['Align', 'Desc', 'Score', 'Start', 'Enzyme', 'Mod', 'Mod2', 'ID'])
i = 0
# df_enz = pd.DataFrame(columns=["Enz", "Seq"])

if sequence=="":
    sys.exit()
else:
    try:
        modstatus = st.status("Finding consensus sites...", state="running", expanded=True)
        for seq in seqs:
          enzyme = df_mods.loc[df_mods['Seq'] == seq, 'Enzyme'].iloc[0]
          mod = df_mods.loc[df_mods['Seq'] == seq, 'Mod'].iloc[0]
          idnum = df_mods.loc[df_mods['Seq'] == seq, 'ID'].iloc[0]
          for a in aligner(refseq, seq, 2, -1, -4, -1):
            # df_enz.loc[len(df_enz.index)] = [enzyme, seq]
            align_start = a[3]
            align_score = a[2]
            i = i+1
            if int(format_alignment2(*a)) == 0:
                pass
            else:
                mod_site = int(format_alignment2(*a))
                accuracy_score = round(align_score*100/(len(seq)*2), 2)
                try:
                    if (refseq[mod_site-1] == seq[20] and accuracy_score >= 60):
                #      desc = f" rRNA modification site is detected with an accuracy percentage of {accuracy_score}%\n rRNA modification base is: {refseq[align_start+20]}({align_start+21})\n"
                      desc = accuracy_score
                      new_row = pd.DataFrame({'Align':[a],'Desc':desc, 'Score':accuracy_score, 'Start':mod_site, 'ModSite':mod_site, 'Enzyme':enzyme, 'Mod':mod, 'Mod2':mod, 'ID':idnum}, index=[0])
                      high_conf_df = pd.concat([new_row,high_conf_df.loc[:]]).reset_index(drop=True)
                
                    elif (refseq[align_start+20] == seq[20] and (5 <= accuracy_score < 60)):
                #      desc = f" rRNA modification site is detected with an accuracy percentage of {accuracy_score}%\n rRNA modification base is: {refseq[align_start+20]}({align_start+21})\n"
                      desc = accuracy_score
                      new_row2 = pd.DataFrame({'Align':[a],'Desc':desc, 'Score':accuracy_score, 'Start':mod_site, 'ModSite':mod_site, 'Enzyme':enzyme, 'Mod':mod, 'Mod2':mod, 'ID':idnum}, index=[0])
                      low_conf_df = pd.concat([new_row2,low_conf_df.loc[:]]).reset_index(drop=True)
                except IndexError:
                    pass
    except IndexError:
        modstatus.update(label="No potential rRNA modification sites were found", state="error", expanded=False)
        st.warning("No alignments could be found that suggests the presence of rRNA modification sites.")
        sys.exit()


#st.write(high_conf_df[['Start','Score', 'Mod']])
# In[20]:
#print(high_conf_df[['Start', 'Enzyme', 'Mod']])


high_conf_df = high_conf_df.sort_values(by=['Score'], ascending=False)
high_conf_df = high_conf_df.groupby(['Start', 'Mod2']).first()
high_conf_df = high_conf_df.sort_values(by=['Score'], ascending=False)

low_conf_df = low_conf_df.sort_values(by=['Score'], ascending=False)
low_conf_df = low_conf_df.groupby(['Start', 'Mod2']).first()
low_conf_df = low_conf_df.sort_values(by=['Score'], ascending=False)


#st.write(low_conf_df[['Score', 'Mod']])
# In[21]:
#print(high_conf_df[['Enzyme', 'Mod']])





# In[23]:


#from contextlib import contextmanager, redirect_stdout
#from io import StringIO
#
##sys.stdout.write("033[1;31m")
##
#@contextmanager
#def st_capture(output_func):
#    with StringIO() as stdout, redirect_stdout(stdout):
#        old_write = stdout.write
#        def new_write(string):
#            ret = old_write(string)
#            output_func(stdout.getvalue())
#            return ret
#        
#        stdout.write = new_write
#        yield
#
#
#output = st.empty()
#with st_capture(output.code):
#    print("High Confidence Hits\n*************************\n")
#    for i in high_conf_df.index:
#        print(format_alignment2(*(high_conf_df['Align'][i])))
#        print(high_conf_df['Desc'][i])
#        print(high_conf_df['Enzyme'][i])
#        print("____________________________________________________________________________\n")
#
#    print("Low Confidence Hits\n*************************\n")
#    for i in low_conf_df.index:
#        print(format_alignment2(*(low_conf_df['Align'][i])))
#        print(low_conf_df['Desc'][i])
#        print(low_conf_df['Enzyme'][i])
#        print("____________________________________________________________________________\n")



# In[23]:
# Color tags
import seaborn as sns
snscols = sns.color_palette("pastel", n_colors=25)
pastels = snscols.as_hex()
mods = df_mods['Mod'].unique()

colors_mod = dict(zip(mods, pastels))

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, ExactPosition, SimpleLocation
import base64
import re

# Create a sequence
sequence_string = sequence
sequence_object = Seq(sequence_string)
 
# Create a record
record = SeqRecord(sequence_object,
                   id='123456789', # random accession number
                   name='Example',
                   description='An example GenBank file generated by BioPython',
                   annotations={"molecule_type": "DNA"})

record_gb = SeqRecord(sequence_object,
                   id=f'{sequence[:5]}', # random accession number
                   name='rRNA_Mods',
                   description='A GenBank file generated in ModrRNA1',
                   annotations={"molecule_type": "DNA"})
 
# Add annotation


if len(high_conf_df) > 0:
    modstatus.update(label="Plotting the sequence chart...", state="running", expanded=True)
else:
    modstatus.update(label="No possible rRNA modification sites were found", state="error", expanded=False)
    sys.exit()

high_conf_df_drop = high_conf_df.drop(high_conf_df.columns, axis=1)

for i in high_conf_df.index:
    mod_point = int(high_conf_df['ModSite'][i])
    mod_type = high_conf_df['Mod'][i]
    mod_id = high_conf_df['ID'][i]
    ac_score = high_conf_df['Score'][i]
    accu = high_conf_df['Desc'][i]
    high_conf_df_drop.loc[i, 'Modification Point'] = mod_point
    high_conf_df_drop.loc[i, 'Modification Type'] = mod_type
    high_conf_df_drop.loc[i, 'Accuracy Score'] = ac_score
    feature = SeqFeature(FeatureLocation(mod_point, mod_point),
                         type = mod_type, 
                         qualifiers = {"label": (f"{mod_type}\n({accu}%)\n{sequence[mod_point-1]}({mod_point})\n{mod_id}")})
    feature2 = SeqFeature(FeatureLocation(mod_point, mod_point),
                         type = mod_type, 
                         qualifiers = {"label": (f"{mod_id}")})
    record.features.append(feature)
    record_gb.features.append(feature2)


high_conf_df_drop = high_conf_df_drop.sort_values(by='Modification Point')
#st.write(high_conf_df_drop



#len_seq = len(sequence)
#seq_code = sequence[:4] + sequence[int(len_seq/2):int(len_seq/2)+4] + sequence[len_seq-4:]

seq_code = hashgen.generate_unique_id(sequence, length=9)

def save_gb():
    genbank_filename = f"./gbs/{seq_code}.gbk"
    
    with open(genbank_filename, "w") as output_file:
        SeqIO.write(record_gb, output_file, "genbank")
    
    with open(f'./gbs/{seq_code}.gbk', 'r') as genbank_file:
        genbank_data = genbank_file.read()
    
    modified_genbank_data = re.sub(r'(\d+)\^(\d+)', r'\1', genbank_data)
    
    
    b64 = base64.b64encode(modified_genbank_data.encode()).decode()
    href = f"data:file/genbank;base64,{b64}"
    down_gb = f'<a href="{href}" download="results.gb">GenBank</a>'
    return down_gb



# In[23]:
from dna_features_viewer import BiopythonTranslator

class CustomTranslator(BiopythonTranslator):

    # Label fields indicates the order in which annotations fields are
    # considered to determine the feature's label
    label_fields = ["label", "note", "name", "gene"]

    def compute_feature_legend_text(self, feature):
        return feature.type
    
    def compute_feature_color(self, feature):
        return colors_mod[feature.type]
    
    def compute_feature_box_color(self, feature):
        return "white"
    
    def compute_feature_box_linewidth(self, feature):
        return 0

def features_properties(f):
    """Mutations get a red label, other features get a pastel color."""
    label = f.qualifiers["label"]
    boxcol = colors_mod[f.type]
#    if f.type == "m5C":
        
#        box_color = "#ffd383"
    return dict(label=label, box_color=boxcol)


 # TRANSLATE AND PLOT THE RECORD

translator = CustomTranslator(features_properties=features_properties)


def plotMods(record):
    graphic_record = translator.translate_record(record)
    fig, axes = graphic_record.plot_on_multiple_lines(
        nucl_per_line=90, plot_sequence=True
    )
    graphic_record.plot_legend(ax=fig, loc=1, ncol=3, frameon=False)
    return fig
    

col1, col2, col3 = st.columns(3)

with col1:
    st.write(' ')

with col2:
    st.image('https://i.ibb.co/2Pyqb37/modlegend.png')

with col3:
    st.write(' ')
    
with st.expander(expanded = True, label = "Potential RNA modification sites:"):
    st.pyplot(plotMods(record))
    if len(high_conf_df_drop) != 0:
        csv_data2 = high_conf_df_drop.to_csv(index=False)
        b64 = base64.b64encode(csv_data2.encode()).decode()
        href = f"data:file/csv;base64,{b64}"
        down_csv = f'<a href="{href}" download="results.csv">CSV</a>'
        down_gb = save_gb()
        st.markdown(f"Download {down_csv} | {down_gb}", unsafe_allow_html=True)



modstatus.update(label="Modification sites found", state="complete", expanded=False)


st.write("""
***
""")

st.subheader('Enzyme Explorer Beta')

taxbact = []
modsbact = []

with open("taxnames", 'r') as file:
    lines = file.readlines()

taxbact = [line.strip() for line in lines]

for i in high_conf_df.index:
    mod_point = int(high_conf_df['ModSite'][i])
    mod_type = high_conf_df['Mod'][i]
    mod_id = high_conf_df['ID'][i]
    modsbact.append(f"{mod_id} - {mod_type}({mod_point})")
    
if "hits" not in st.session_state:
    st.session_state["hits"] = None

def hitsUp():
    st.session_state["hits"] = None


with st.form(key='my_form2'):
    option1 = st.selectbox(
       "Select the source bacteria of this sequence: (if known)",
       (taxbact),
       index=None,
       placeholder="Search bacteria",
    )
    option2 = st.selectbox(
       "Select the modification site code",
       (modsbact),
       index=None,
       placeholder="Search the modification",
    )
    submit_button2 = st.form_submit_button(label='Submit')
    if submit_button2: print("submit")
    if option1 and option2:
        site = option2.split(" - ")[0] if " - " in option2 else option2
        caption = option2.split(" - ")[1] if " - " in option2 else option2
        enz = high_conf_df.loc[high_conf_df['ID'] == site, 'Enzyme'].iloc[0]
    else:
        sys.exit()
    



#st.code(high_conf_df['Enzyme'])


if option1 and option2:
    if (enz == None):
        st.error("Enzyme Explorer is not available for this site.")
        sys.exit()
    else:
        st.write('You selected:', option1)
        


from ncbi import ete3
ncbi = ete3.ncbiquery.NCBITaxa()

if option1 and option2:
    taxdict = ncbi.get_name_translator([option1])
    
    first_key = next(iter(taxdict))
    first_value = taxdict[first_key]
    taxid = int(first_value[0])
    
    #print(taxid)
    
    st.write('TaxID:', taxid)

# In[9]:

def first_lowest(df):
    min_value = df.min()
    row_index = (df == min_value).idxmax(axis=0)
    return row_index

def multiplyList(myList):
    result = 1
    for x in myList:
        result = result * x
    return result

from requests import get, post
from time import sleep

enz_pdb = f"pdbs/{enz}.pdb"


with open(enz_pdb,'r') as file:
    pdbs = file.read()



modes = ["pdb100","afdb50","afdb-swissprot","afdb-proteome"]
tickets = []

@st.cache_data
def submitFold(site):
    print(site)
    for i in range(len(modes)):
        tickets.append(post('https://search.foldseek.com/api/ticket', {
                'q' : pdbs,
                'database[]' : modes[i],
                'mode' : '3diaa',
            }).json())
    # poll until the job was successful or failed
    repeat = True
    status = []
    complete=0
    try:
        while repeat:
            for t in tickets:
                apnd = get('https://search.foldseek.com/api/ticket/' + t['id']).json()
                if apnd not in status:
                    status.append(apnd)
            for s in status:
                if s['status'] == "ERROR":
                # handle error
                    st.error("Error when searching!")
                    sys.exit(0)  
        
            print(status)
            # wait a short time between poll requests
            sleep(1)
            for w in status:
                if w['status'] == "COMPLETE":
                    complete+=1                
            repeat = (complete<4)       
    except Exception:
        # This block will be executed when any exception occurs
        st.warning('Please try hitting Submit again.')
    return tickets 
    
# get all hits for the first query (0)
hits=[]


@st.cache_data
def getFold(t):
    return get('https://search.foldseek.com/api/result/' + t['id'] + '/0').json()


if option1 and option2:
    try:
        with st.spinner('Please Wait...'):
            for t in submitFold(site):
                result = getFold(t)
                hits.append(result)       
    except Exception:
        # This block will be executed when any exception occurs
        st.warning('Please try hitting Submit again.')

# # print pairwise alignment of first hit of first database
# print(result['results'][0]['alignments'][0]['qAln'])
# print(result['results'][0]['alignments'][0]['dbAln'])

def getTarget(text):
    parts = text.split("_")
    first_part = parts[0]
    final_parts = first_part.split("-")
    if len(final_parts) > 1:
        result = "-".join(final_parts[:2])
    else:
        result = final_parts[0]
    return result


#st.code(len(st.session_state["hits"]))
if not len(hits) == 0:
    st.session_state["hits"] = hits

#st.code(len(st.session_state["hits"]))

pdb100_hits = hits[0]['results'][0]['alignments']
afdb50_hits = hits[1]['results'][0]['alignments']
swissprot_hits = hits[2]['results'][0]['alignments']
proteome_hits = hits[3]['results'][0]['alignments']

combined_hits = pdb100_hits + afdb50_hits + swissprot_hits + proteome_hits


enz_df = pd.DataFrame(columns=['protein','distance', 'prob', 'qscore', 'align_seq'])
dist_ecoli = ncbi.get_topology([str(taxid), '83333'], intermediate_nodes=True)


#st.code(pdb100_hits[1].keys())
#st.code(afdb50_hits[2]['tSeq'])
#align_seq = (afdb50_hits[3]['tSeq'])[afdb50_hits[3]['dbStartPos']-1:afdb50_hits[3]['dbEndPos']]
#st.code(align_seq)

if "folded" not in st.session_state:
    st.session_state["folded"] = None

if "enz_df" not in st.session_state:
    st.session_state["enz_df"] = None


@st.cache_data
def FoldSeek(hits, taxid):
    st.session_state["folded"] = True
    for i in range(len(combined_hits)):
    #    if combined_hits[i]['taxId'] == taxid:
    #        st.write(f"Putative enzymes for {caption} modification of {combined_hits[i]['taxName']} rRNA are:")
    #        target = combined_hits[i]['target'][0:4]
    #        st.markdown(f"https://www.rcsb.org/pdb/explore.do?structureId={target}")
    #        break
    #    else:
            target = combined_hits[i]['target']
            prob = combined_hits[i]['score']
            qalign = combined_hits[i]['qEndPos'] - combined_hits[i]['qStartPos'] 
            qlen = combined_hits[i]['qLen']
            qscore = round((qalign*100)/qlen, 0)
            align_seq = (combined_hits[i]['tSeq'])[combined_hits[i]['dbStartPos']-1:combined_hits[i]['dbEndPos']]
#            if "tRNA" in target:
#                break
            if combined_hits[i]['taxId'] == 0 or combined_hits[i]['taxId'] == 1094892:
                distance = ncbi.get_topology([str(taxid), '9606'], intermediate_nodes=True)
                enz_df.loc[len(enz_df)] = {'protein':target,'distance':distance, 'prob': prob, 'qscore' : qscore, 'align_seq':align_seq}
            else:
                distance = ncbi.get_topology([str(taxid), str(combined_hits[i]['taxId'])], intermediate_nodes=True)
                enz_df.loc[len(enz_df)] = {'protein':target,'distance':distance, 'prob': prob, 'qscore' : qscore, 'align_seq':align_seq}
    return enz_df.loc[(enz_df['qscore'] >= 50) & (enz_df['prob'] >= 60)]  

  
#enz_df = enz_df[enz_df['distance'] < dist_ecoli]

#enz_hits = FoldSeek(combined_hits)

try:
    with st.spinner('Please Wait...'):                
        st.session_state["enz_df"] = FoldSeek(combined_hits, taxid)
except Exception:
    # This block will be executed when any exception occurs
    st.warning('Please try hitting Submit again.')
#st.toast("FoldSeek Done!")

st.session_state["enz_df"].drop(st.session_state["enz_df"][st.session_state["enz_df"]['protein'].str.contains('tRNA')].index, inplace=True)    
            




left, mid1, mid2, right = st.columns(4, gap="large")
with left: 
    on = st.toggle('Sort by Enzyme Similarity')
    if on:
        result_enz = st.session_state["enz_df"].sort_values(by=['prob','distance'], ascending = [False, True])
    else:
        result_enz = st.session_state["enz_df"].sort_values(by=['distance','prob'], ascending = [True, False])
with right:
    limit = st.select_slider('Hits limit:', options = [1,5,10,15,20,25], value = 5, label_visibility="collapsed")

result_enz = result_enz.reset_index(drop=True)[:limit]

result_enz['link'] = None

for i in range(len(result_enz)):
    target_cut = getTarget(result_enz['protein'].iloc[i])
    if len(target_cut) == 4:
        link = f"https://www.rcsb.org/pdb/explore.do?structureId={target_cut}"
    else:
        link = f"https://www.alphafold.ebi.ac.uk/entry/{target_cut}"
    result_enz.loc[i, 'link'] = link

result_enz_no = result_enz

def make_clickable(link):
    return f'<a target="_blank" href="{link}">Go to Entry</a>'

# link is the column with hyperlinks
result_enz['link'] = result_enz['link'].apply(make_clickable)
result_enz2 = result_enz.drop('align_seq', axis=1).to_html(escape=False, justify="center")

st.caption("If top hits does not seem sensible (not containing keywords like rRNA and the respective modification type), try sorting by enzyme similarity.")
with st.expander("Sorted in order of most related to least related:", expanded=True):
#    st.markdown(result_enz2, unsafe_allow_html=True)    
    st.markdown(
        f"""
        <div style="display: flex; justify-content: center;">
            <div>
                {result_enz2}
            </div>
        </div>
        """,
        unsafe_allow_html=True,
    )

def convert_df(df):
    return df.to_csv().encode('utf-8')


count = 0
allhits = ""
for t in submitFold(site):
    allhits = allhits + f"[{modes[count]}](https://search.foldseek.com/result/{t['id']}/0) | "
    count = count + 1
#if "option" not in st.session_state:
#    st.session_state["align_seq"] = None



if not len(result_enz) == 0:
    st.caption("Less distance, higher prob, better match.")
    csv_data = result_enz_no.to_csv(index=False)
    b64 = base64.b64encode(csv_data.encode()).decode()
    href = f"data:file/csv;base64,{b64}"
    down_csv = f'<a href="{href}" download="mydata.csv">:green[Download above results as CSV]</a>'
    st.markdown(f"All hits: {allhits}\t\t{down_csv}", unsafe_allow_html=True)
    

st.write("""
***
""")

st.subheader('Protein BLAST')

st.markdown(f"You can BLAST any of these enzymes on ***{option1}*** genome(s). :red[(This may take a while.)]")


with st.form(key='my_form3'):
    option3 = st.selectbox(
       "Select the enzyme",
       (result_enz['protein'].tolist()),
       index=0,
       placeholder="Search the enzyme",
    )
    left, col1, col2, col3, col4 = st.columns(5, gap="large")
    with left:
        option4 = st.selectbox(
       "Which BLAST?",
       (['blastp', 'tblastn']),
       index=0,
       placeholder="Search the enzyme",
    )
    submit_button3 = st.form_submit_button(label='BLAST!')    
    if option3 and option4:
        enz_seq = result_enz.loc[result_enz['protein'] == option3, 'align_seq'].iloc[0]
        blastop = option4
    else:
        sys.exit()
#    enz = high_conf_df.loc[high_conf_df['ID'] == site, 'Enzyme']
    


if enz_seq is not None and submit_button3:
    st.caption("If blast() is taking too long, please be kind enough use the below sequence to do a BLAST yourself. FoldSeek aligned sequence:")
    st.code(enz_seq)


from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

pdb_seq = enz_seq

print("done")

@st.cache_data    
def blast(pdb_seq, taxid, blastop):
    if blastop == 'tblastn':
        return NCBIWWW.qblast('tblastn', 'nt', pdb_seq, entrez_query=f"txid{str(taxid)}[ORGN]")
    else:
        return NCBIWWW.qblast('blastp', 'nr', pdb_seq, entrez_query=f"txid{str(taxid)}[ORGN]")

def delayed_print():
    st.warning("If blast() is taking too much time, please use the above sequence to do a BLAST yourself.")


if enz_seq is not None and submit_button3:
    result_handle = blast(pdb_seq, taxid, blastop)
else:
    sys.exit()
    
# Parse and print BLAST results
blast_records = NCBIXML.parse(result_handle)


hit = 0
for blast_record in blast_records:
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            hit = hit+1
            with st.expander(label=f"Hit {hit}", expanded=True):
                st.write(f"{alignment.title}")
                st.code(f"""{hsp}""")
                accession = alignment.accession
                if blastop == 'blastp':
                    st.write(f"<a href='https://www.ncbi.nlm.nih.gov/protein/{accession}' id='my-link'>{accession}</a>", unsafe_allow_html=True)
                elif blastop == 'tblastn':
                    st.write(f"<a href='https://www.ncbi.nlm.nih.gov/nuccore/{accession}' id='my-link'>{accession}</a>", unsafe_allow_html=True)
                    

if hit == 0:
    st.warning(f"No matches found in {option1} genome")








# In[42]:



#

#

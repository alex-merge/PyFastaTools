# -*- coding: utf-8 -*-
"""
Fasta file tools for python

Created on Fri Mar 18 08:26:12 2022

@author: Alex-932
@version: 1.0
"""

import pandas as pd
import numpy as np
import time
import os
import re

class PyFastaTools() :
    
    def __init__(self, path=os.getcwd()):
        """
        Initialize the dataframe from the given file.

        Parameters
        ----------
        path : str
            Path to the working folder.

        Returns
        -------
        None.

        """
        if os.path.isfile(path) :
            if PyFastaTools.is_fasta(path):
                self.file_path = path
            else :
                raise TypeError("The selected file is not supported !")
        else :
            self.file_path = PyFastaTools.file_explorer(path)
        self.df = PyFastaTools.df_creator(self.file_path)
        
    def file_explorer(path):
        """
        Search in an iterative way for a file and ask the user to choose 
        a file.

        Parameters
        ----------
        path : str
            Path to search in.

        Returns
        -------
        str
            Path to the selected file.

        """
        present_file = os.listdir(path)
        
        for k in range(len(present_file)):
            print(k, " - ", present_file[k])
        
        choice = int(input("What file (index) : "))
        
        new_path = path+'\\'+present_file[choice]
        if os.path.isdir(new_path) :
            return PyFastaTools.file_explorer(new_path)
        elif not PyFastaTools.is_fasta(new_path) :
            raise TypeError("The selected file is not supported !")
            return PyFastaTools.file_explorer(path)
        return new_path
        
    def is_fasta(path):
        splitted_path = re.split(r'\\',path)
        extension = re.split(r'\.',splitted_path[-1])
        if extension[-1] in ["fq","fastq","fasta","fa"] :
            return True
        return False
    
    def checker(string, allowed_char):
        """
        Method that check if the given string has unallowed characters in it.
    
        Parameters
        ----------
        string : str
            The string that will be checked.
        allowed_char : list
            List taht contains the allowed characters.
    
        Returns
        -------
        bool
            False if the string contains an unallowed character.
            True if not.
    
        """
        check = list(string)
        for k in check :
            if k not in allowed_char :
                return False
        return True
    
    def cleanup(file_list, element='\n' , position=0):
        """
        Remove the "/n" that are in the end of all strings within the list 
        coming from the .readlines() method.
    
        Parameters
        ----------
        file_list : list
            List that come from the .readlines() method.
    
        Returns
        -------
        list
            List of strings without the "/n" at the end of them.
    
        """
        return [re.split(element, k)[position] for k in file_list]
                    
    def df_creator(path):
        """
        Create a pandas dataframe with columns representing : header, sequence 
        and the quality in the case of a FastaQ file.
    
        Parameters
        ----------
        file : str
            Path to the file.
        element : str, optional
            Element that will split the string. Default is '\n'.
        position : int, optional
            Position of the searched sequence in the list given by the method
            re.split(). Default is 0.
    
        Returns
        -------
        df : pandas.dataframe
            Dataframe containing the informations of the file.
    
        """
        fastafile = open(path)
        line_list = PyFastaTools.cleanup(fastafile.readlines()) 
        #list that contains all the lines of the given file.
        
        #We're going to create a list of all "paragraphs"
        #We first need to recognize the index of the sequence line
        seq_index = [k for k in range(len(line_list)) \
                        if PyFastaTools.checker(line_list[k],\
                                                ['A','T','G','C','N'])]
        seq_index.append(len(line_list)+1)
        
        #Now we add the paragraphs into one list
        par_list_raw = [line_list[seq_index[k]-1:seq_index[k+1]-1] \
                    for k in range(len(seq_index)-1)]
            
        par_list = []
        for par in par_list_raw:
            if len(par) == 4:
                #If the paragraph is 4 lines long then it's a FastaQ so the 3rd
                #line is useless as it is a recall of the header
                par[0] = PyFastaTools.cleanup([par[0]], element=' ', \
                                              position=0)[0]
                par.pop(2)
                par_list.append(par)
            else :
                #Else, it's a Fasta file so no line deleting but we add a 
                #quality score of np.nan
                par.append(None)
                par_list.append(par)
        
        #We create the pandas dataframe 
        df_labels = ["Header","Sequence","Quality"]
        df = pd.DataFrame(np.array(par_list),\
                          index=[k for k in range(1,len(par_list)+1)],\
                              columns=df_labels)
        
        #Then we add some informations concerning the reads to the dataframe 
        seq_series = df["Sequence"]
        
        #In the case of a fasta file there is no quality score so we put None
        #in place of the actual quality score per base list.
        tempory_list = []
        for k in df["Quality"] :
            if k != None :
                tempory_list.append(PyFastaTools.to_Qscore(k))
            else :
                tempory_list.append(None)
        df["QScore"] = tempory_list
        
        #Also in the case of a None type quality, we do not calculate the mean.
        #But we set it to None.
        tempory_list = []
        for k in df["QScore"] :
            if k != None :
                tempory_list.append(np.mean(k))
            else :
                tempory_list.append(None)
        df["QScore_Mean"] = tempory_list
        
        #Then we add some more informations regarding the sequence.
        df["Length"] = [len(k) for k in seq_series]
        df["N_count"] = [k.count("N") for k in seq_series]
        df["GC_proportion"] = [(k.count("G")+k.count("C"))/len(k) \
                                for k in seq_series]
            
        return df
    
    def fasta_summary(self):
        """
        Give a summary of the given dataframe.
    
        Returns
        -------
        None.
    
        """
        print("Number of reads : ", self.df.shape[0])
        print("Average Q-scores : ", self.df.QScore_Mean.mean())
        print("Average length of the sequences : ", self.df.Length.mean())
        print("Median of the sequence length : ", self.df.Length.median())
        print("Average N in the sequences : ", self.df.N_count.mean())
        print("Average proportion of G and C in the sequences : ",\
              self.df.GC_proportion.mean())
            
    def search_header(self, header):
        """
        Search for the sequence of a given header.
    
        Parameters
        ----------
        dataframe : pandas.DataFrame
            Dataframe to search into.
        header : str
            Header to search for.
    
        Returns
        -------
        pandas.Series
            Series with the sequence that correspond to the given header.
    
        """
        return self.df.loc[ self.df["Header"] == header ]\
            ["Sequence"]
        
    def search_length(self, length, mode='='):
        """
        Search for rows where the length of the sequence respect the given
        conditions.
    
        Parameters
        ----------
        dataframe : pandas.DataFrame
            Dataframe to search into.
        length : int
            Length condition.
        mode : str, optional
            Search condition between '=', '>', '<', '!='. The default is '='.
    
        Returns
        -------
        pandas.DataFrame
            Dataframe with the header and the sequence of the rows
            that did respect the given conditions.
    
        """
        if mode == '=' :
            return self.df.loc[ self.df["Length"] == length ]\
                [["Header","Sequence"]]
        if mode == '>' :
            return self.df.loc[ self.df["Length"] > length ]\
                [["Header","Sequence"]]
        if mode == '<' :
            return self.df.loc[ self.df["Length"] < length ]\
                [["Header","Sequence"]]
        if mode == '!=' :
            return self.df.loc[ self.df["Length"] != length ]\
                [["Header","Sequence"]]
    
    def search_pattern(self, pattern):
        """
        Use the pandas.str.contains() method to search for rows that contains 
        the given pattern within its sequence.
    
        Parameters
        ----------
        dataframe : pandas.DataFrame
            Dataframe to search into.
        pattern : str
            Pattern to search for.
    
        Returns
        -------
        pandas.DataFrame
            Dataframe that contains the rows where there has been a match.
    
        """
        return self.df.loc[ self.df["Sequence"].str.contains(pattern) ]
    
    def to_Qscore(string):
        """
        Convert the Quality string from a symbol to the Q-Score.
    
        Parameters
        ----------
        string : str
            Quality with symbols in one string.
    
        Returns
        -------
        list
            List that contains the Q-Score of each base of a read.
    
        """
        #The ord() method directly provide the ASCII value of a symbol
        return [ord(k)-33 for k in list(string)]
    
    def quality_check(self, min_length=80, av_qscore_min=30, min_qscore=30,\
            adapterSeqList=['AGATCGGAAGAGC','TGGAATTCTCGG','CTGTCTCTTATA']):
        """
        Method to clean a fasta file based on quality. 
    
        Parameters
        ----------
        min_length : Int, optional
            Minimal length that reads must have to pass. The default is 80.
        av_qscore_min : Int, optional
            Average Q-Score that reads must have to pass. The default is 30.
        min_qscore : Int, optional
            Q-Score that a base must have to pass. The default is 30.
    
        Returns
        -------
        df_qc : pandas.dataframe
            Return a cleaned dataframe.
    
        """
        
        df_qc = self.df.copy()
        #Every iteration, rows that contains the screened adapter are removed.
        for k in adapterSeqList :
            df_qc = df_qc.loc[self.df["Sequence"].str.contains(k) == 0]
        
        #Every base that do not pass the check is replaced with a "N".
        clean_reads = []
        validation_list = df_qc["QScore"].isnull()
        for index in range(1,len(df_qc["QScore"])+1) :
            #"index" represent the row's index in the dataframe.
            if not validation_list[index] :
                cr = []
                for value in range(len(df_qc["QScore"][index])):
                    #"value" represent the index in the Q-Score list.
                    #The aim here is to rebuild the sequence base per base.
                    if df_qc["QScore"][index][value] >= min_qscore :
                        #If it pass the check, we add the base to the 
                        #new sequence. 
                        cr.append(list(df_qc["Sequence"][index])[value])
                    else :
                        #If it doesn't, the base is replaced by a "N" in the 
                        #new seq.
                        cr.append("N")
                clean_reads.append("".join(cr))
            else : 
                clean_reads.append(df_qc["Sequence"][index])
            
        #We replace the sequences with the corrected ones in the case of a 
        #FastQ paragraph.
        df_qc["Sequence"] = clean_reads
        
        #First the reads that are long enough are kept.
        df_qc = df_qc.loc[df_qc["Length"] > min_length]
        #Then those with an average Q-Score greater than the limit are kept. 
        df_qc = df_qc.loc[df_qc["QScore_Mean"] > av_qscore_min]
        
        #Finally, the "N" within the sequences are reprocessed.
        df_qc["N_count"] = [k.count("N") for k in df_qc["Sequence"]]
        
        self.df_qc = df_qc
        
        return self.df_qc
            
    
    
if __name__ == "__main__":
    df = PyFastaTools()
    start = time.time()
    df.fasta_summary()
    df.quality_check()
    print( "Execution time : {} s".format(time.time()-start) )

    
    

    

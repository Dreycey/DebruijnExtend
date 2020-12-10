import sys
import pickle
import matplotlib.pyplot as plt
import numpy as np 

def plot_confusionmatrix():                                                        
    """                                                                         
    This function plot the 
    """ 

def avg_acc_per_length(dict_list, i):
    """
    This plots a scatter plot for a particular software. The scatterplot 
    shows accuracy vs length
    """
    N = 3
    total_accuracy_array = []
    total_length_array = []
    software = dict_list[i]
    # combine data
    for k_fold in range(N):
        accuracy_array = software[k_fold][0]
        length_array = software[k_fold][1]
        total_accuracy_array.append(accuracy_array)
        total_length_array.append(length_array)
    # plot
    plt.scatter(total_length_array, total_accuracy_array, color="orange", edgecolor="b")
    plt.ylim(0, 1)
    plt.ylabel("Accuracy")
    plt.xlabel("Protein Length")
    plt.show()

def difference_hist(dict_list, i, j):
    """
    This function plots a histogram of the difference in accuracy between two
    competing softwares i and j.
    """
    N = 3
    software_one = dict_list[i]
    software_two = dict_list[j]
    # calculate distance
    difference_array = []
    for k_fold in range(N):
        difference_array.append([])
        for protein in range(len(software_one[k_fold][0])):
            difference = abs(software_one[k_fold][0][protein] - 
                             software_two[k_fold][0][protein])
            difference_array[k_fold].append(difference)

    # plot the difference 
    for k_fold in range(N):
        plt.hist(difference_array[k_fold], color="orange", edgecolor="black")
        title = "Differences in Protein to Protein Accuracy"
        tool_comp = "DebruijnExtend vs PsiPred"
        plt.title(f"{title} \n {tool_comp} \n Test Set {k_fold+1}")
        plt.ylabel("# of proteins")
        plt.xlabel("Difference in Accuracy")
        plt.show()

def plot_bar(dict_list):
    """
    This method plots a simple bar graph for the data.
    input: list of dicts per software.
    """
    N = 3
    fig, ax = plt.subplots()
    # create array with means and STDs
    meadSTD_dict = {} 
    for software_i in range(len(dict_list)):
        meadSTD_dict[software_i] = [[],[]] # for mean and std
        for k_fold in range(N):
            accuracy_array = dict_list[software_i][k_fold][0]    
            mean = np.mean(accuracy_array)
            std = np.std(accuracy_array)
            meadSTD_dict[software_i][0].append(mean)
            meadSTD_dict[software_i][1].append(std)
    
    # create bar plots
    width = 0.35
    ind = np.arange(N) 
    bar_array = []
    for software_i in range(len(dict_list)):
        p_temp = ax.bar(ind + width*software_i, 
                        meadSTD_dict[software_i][0], 
                        width, 
                        bottom=0, 
                        yerr=meadSTD_dict[software_i][1])
        bar_array.append(p_temp)
    bar_tuple = tuple(bar_array)

    # plot the bar plot
    ax.set_title('Secondary Structure Prediction Accuracy')
    ax.set_ylabel('Accuracy')
    ax.set_xticks(ind + width / 2)
    ax.set_xticklabels(('Test Set 1', 'Test Set 2', 'Test Set 3'))

    ax.legend((bar_array[0][0], bar_array[1][0]), ('DebruijnExtend', 'PsiPred'))
#    ax.yaxis.set_units(inch)
#    ax.autoscale_view()
    plt.show()

def main():
    # add each software dict to a list.
    dict_list = []
    for soft_i in range(1,len(sys.argv)):
        dictionary = open(sys.argv[soft_i], "rb")
        opened_pickle = pickle.load(dictionary)
        dict_list.append(opened_pickle)
        dictionary.close() # close the file opened.

    # create a bar plot.
    plot_bar(dict_list)
    
    # plot the difference histogram
    difference_hist(dict_list, 0, 1)

    # accuracy vs length
    avg_acc_per_length(dict_list, 0)
    avg_acc_per_length(dict_list, 1)
if __name__ == "__main__":
    main()

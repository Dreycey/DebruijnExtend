import sys
import pickle


def plot_average_perk():
    """
    """


def plot_heatmap():                                                        
    """                                                                         
    """ 

def plot_hist():
    """
    """

def avg_acc_per_length():
    """
    """

def rms_difference():
    """
    """

def main():
    dictionary = open(sys.argv[1], "rb")
    opened_pickle = pickle.load(dictionary)
    for k, v in opened_pickle.items():
        print(k, v)
        print("\n")
    opened_pickle.close()
if __name__ == "__main__":
    main()

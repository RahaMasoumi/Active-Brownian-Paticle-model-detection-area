import pandas as pd
import numpy as np
from plotly.subplots import make_subplots
import plotly.graph_objects as go


def conversion(path):
    """
    Converts a tij.dat file into a np.array.

    :param path: path of the tij.dat file
    :type path: str
    :return: np.array of the tij data
    :rtype: np.array
    """
    df = pd.read_csv(path, sep='\t')
    tij_array = df.to_numpy()
    return tij_array


def unique(ar):
    """
    This function gives each unique value of an array and the number of occurrence of the value

    :param ar: Array that is studied
    :type ar: np.array
    :return: Unique values of ar and the number of occurrences
    :rtype: tuple of np.array
    """
    values, counts = np.unique(ar, return_counts=True)
    return values, counts


def common(ar1, ar2):
    """
    This functions returns the common rows of ar1 and ar2
    :param ar1: First array
    :type ar1: np.array
    :param ar2: Second array
    :type ar2: np.array
    :return: array of common rows
    :rtype: np.array
    """
    common_array = np.array([x for x in set(tuple(x) for x in ar1) & set(tuple(x) for x in ar2)])
    return common_array


def lost(ar1, ar2):
    """
    This function finds the rows that are in ar1 but not in ar2. These rows are called the lost rows.
    :param ar1: First array
    :type ar1: np.array
    :param ar2: Second array
    :type ar2: np.array
    :return: array of the lost rows
    :rtype: np.array
    """

    set1 = {tuple(x) for x in ar1}
    set2 = {tuple(x) for x in ar2}
    lost_set = (set1 ^ set2) & set1
    if len(lost_set) != 0:
        lost_array = np.array(list(lost_set))
    else:
        lost_array = np.empty((0, 2), dtype=int)
    return lost_array


def new(ar1, ar2):
    """
    This function finds the rows that are in ar2 but not in ar1. These rows are called the new rows.
    :param ar1: First array
    :type ar1: np.array
    :param ar2: Second array
    :type ar2: np.array
    :return: array of the lost rows
    :rtype: np.array
    """
    set1 = {tuple(x) for x in ar1}
    set2 = {tuple(x) for x in ar2}
    new_set = set2 - set1
    if len(new_set) != 0:
        new_array = np.array(list(new_set))
    else:
        new_array = np.empty((0, 2), dtype=int)
    return new_array


def add_time(time, couples, time_sequence_array): #the function was editted on 2023-01-25
    """
    This function adds
    :param time:
    :param couples:
    :param timeline_array:
    :return:
    """
    for elt in couples:
        i = elt[0]
        j = elt[1]
        if i < j:
            time_sequence_array[i, j].append(time)
        else:
            time_sequence_array[j, i].append(time)
    return time_sequence_array


def time_sequence(tij_array):
    """
    This function returns an array of timelines of interactions between all the particles. A timeline between particle
    i and j has the following form [t1, t2, t3, t4 ...] with all the odd elements the time of the beginning of an
    interaction and all the even elements the time of the end of an interaction. As the interaction between i and j is
    strictly the same as the interaction between j and i the array should be symmetric, with all the diagonal elements
    equal to 0 (no interaction between i and i). In our case the array is strictly upper triangular (no need to keep in
    memory all the elements).

    :param tij_array: Array of the tij elements, that are needed to create the timeline array
    :type tij_array: np.array
    :param dt: Increment of time for each step
    :type dt: float or int
    :return: Array of timelines.
    :rtype: np.array of lists
    """
    time_array, counts = unique(tij_array[:, 0]) #the first column is the time of collision
    ij_array = tij_array[:, 1:] #the interaction parts of array (removes the first column=time)
    ij_array = np.int64(ij_array) # numpy.power evaluates 100 ** 8 correctly for 64-bit integers, 
    i_min = np.min(ij_array) #min ij_array
    i_max = np.max(ij_array) #max ij_array

    ij_array = ij_array - i_min  #we do this for making the labels of the nodes starts from 0!

    time_sequence_size = (i_max - i_min + 1,) * 2 # timeline_array [i_max- i_min +1][i_max- i_min +1]
    time_sequence_array = np.frompyfunc(list, 0, 1)(np.empty(time_sequence_size, dtype=object))
    count = counts[0] # the first element of time_array
    couples = ij_array[0:count] #the first row of ij_array and contact between node 180, and 241 
    old_time = time_array[0]
    time_sequence_array = add_time(old_time, couples, time_sequence_array)  




    for step, time in enumerate(time_array[1:]):

        new_count = count + counts[step + 1]
        new_couples = ij_array[count: new_count, :] #the next rows of the ij_array

        time_sequence_array = add_time(time, new_couples, time_sequence_array)



        count = new_count
        old_time = time
    return time_sequence_array

###################################################################################################
def timeline(time_sequence_array ,dt):
    
    for i in range(len(time_sequence_array)):
        
        for j in range(i+1,len(time_sequence_array) ):
            lst=[]  
            elt1=time_sequence_array[i][j]
            
            if len(elt1)>1:
                t1=elt1[0]
                lst.append(t1)
                
                for ii in range(1,len(elt1)):
                    t2=elt1[ii]

                    if t2-t1==dt :
                        t1=t2

                    else:
                        lst.append(t1+dt)
                        lst.append(t2)
                        t1=t2

                    if ii==len(elt1)-1:
                        lst.append(t2+dt)

            elif len(elt1)==1:
                t1=elt1[0]
                lst.append(t1)
                
                t2=t1+dt
                lst.append(t2)
                
            else:
                pass
            
            time_sequence_array[i][j]=lst[:]

    return time_sequence_array

######################################################################################################

def quantities_calculator(timeline_array, dec=1): #what is timeline_array #raha
    """
    Calculates 4 different quantities - contact time, inter-contact time, number of contacts and weight - that are
    needed to compare and validate different models with real data.

    :param timeline_array: Array of timelines.
    :type timeline_array: np.array of lists
    :param dec: decimals to which we around the quantities. Default is equal to 1
    :type dec: int, optional
    """
    contact_time_array = []
    inter_contact_time_array = []
    number_contact_array = []
    link_weight_array = []
    for elt in timeline_array: #what is this array? timeline_array #raha

        for elt1 in elt:
   
            if len(elt1) % 2 == 1: # if the length of the elt1 is odd then we delete the last element with below consrucion  #raha
                elt1.pop()

            if len(elt1) > 0:
                number_contact_array.append(len(elt1) // 2)  #floor division
                contact_time = [b - a for a, b in tuple(zip(elt1, elt1[1:]))[::2]] #why like this?
                contact_time_array.extend(contact_time)
                link_weight_array.append(sum(contact_time))
                inter_contact_time = [b - a for a, b in tuple(zip(elt1[1:], elt1[2:]))[::2]]
                inter_contact_time_array.extend(inter_contact_time)

    contact_time_array, inter_contact_time_array = np.array(contact_time_array), np.array(inter_contact_time_array)
    number_contact_array, link_weight_array = np.array(number_contact_array, dtype=int), np.array(link_weight_array)
    contact_time_array = np.around(contact_time_array, decimals=dec)
    inter_contact_time_array = np.around(inter_contact_time_array, decimals=dec)
    link_weight_array = np.around(link_weight_array, decimals=dec)

    return contact_time_array, inter_contact_time_array, number_contact_array, link_weight_array


def regroup_data(ar):
    """
    This function regroups the quantities with the same value and calculates the number of occurrence of the value.
    The results are then put in a array where for all i, the first element of row i is value i and the second element
    of row i is its number of occurrences.

    :param ar: Array of all the values, of shape (n, )
    :type ar: np.array
    :return: array of shape (n', 2) of values and counts
    :rtype: np.array
    """
    values, counts = unique(ar)
    return np.concatenate((values.reshape((-1, 1)), counts.reshape((-1, 1))), axis=1)


def representation(quantities, title, scale='linear'):
    """
    Represents 4 different quantities - contact time, inter-contact time, number of contacts and weight - in histograms.

    :param quantities: tuple of the 4 quantities that are represented
    :type quantities: tuple of np.arrays
    :param title: Title of the figure
    :type title: str
    :param scale: Scale of the plot. Can be 'linear' (default), 'log' or 'semi-log'
    :type scale: str, optional
    """
    fig = make_subplots(rows=2, cols=2)
    index = [[1, 1], [1, 2], [2, 1], [2, 2]]

    if scale == 'log':
        scale_x, scale_y = scale, scale

    elif scale == 'linear':
        scale_x, scale_y = scale, scale

    else:
        scale_x, scale_y = 'linear', 'log'

    # Update xaxis properties
    fig.update_xaxes(title_text="Contact duration", type=scale_x, row=1, col=1)
    fig.update_xaxes(title_text="Intercontact duration", type=scale_x, row=1, col=2)
    fig.update_xaxes(title_text="Number of contacts", type=scale_x, row=2, col=1)
    fig.update_xaxes(title_text="weight", type=scale_x, row=2, col=2)

    # Update yaxis properties
    fig.update_yaxes(title_text="Distribution of contact duration", type=scale_y, row=1, col=1)
    fig.update_yaxes(title_text="Distribution of intercontact duration", type=scale_y, row=1, col=2)
    fig.update_yaxes(title_text="Distribution of number of contacts", type=scale_y, row=2, col=1)
    fig.update_yaxes(title_text="Distribution of weight", type=scale_y, row=2, col=2)

    for i, data in enumerate(quantities):
        a = index[i][0]
        b = index[i][1]
        if scale == 'log':
            counts, bins = np.histogram(data, bins=np.logspace(np.log10(np.min(data - 0.5)),
                                                               np.log10(np.max(data + 0.5))), density=True)

        else:
            counts, bins = np.histogram(data, bins='auto', density=True)
        bins = 0.5 * (bins[:-1] + bins[1:])
        fig.add_trace(go.Scatter(x=bins, y=counts, mode='markers', showlegend=False), row=a, col=b)

    fig.show()


def make_hist(quantities, title, scale='linear'):
    """
    Represents 4 different quantities - contact time, inter-contact time, number of contacts and weight - in histograms.

    :param quantities: tuple of the 4 quantities that are represented
    :type quantities: tuple of np.arrays
    :param title: Title of the figure
    :type title: str
    :param scale: Scale of the plot. Can be 'linear' (default), 'log' or 'semi-log'
    :type scale: str, optional
    """
    fig = make_subplots(rows=2, cols=2)
    index = [[1, 1], [1, 2], [2, 1], [2, 2]]

    if scale == 'log':
        scale_x, scale_y = scale, scale

    elif scale == 'linear':
        scale_x, scale_y = scale, scale

    else:
        scale_x, scale_y = 'linear', 'log'

    # Update x axis properties
    fig.update_xaxes(title_text="Contact duration", type=scale_x, row=1, col=1)
    fig.update_xaxes(title_text="Inter contact duration", type=scale_x, row=1, col=2)
    fig.update_xaxes(title_text="Number of contacts", type=scale_x, row=2, col=1)
    fig.update_xaxes(title_text="weight", type=scale_x, row=2, col=2)

    # Update y axis properties
    fig.update_yaxes(title_text="Contact duration distribution", type=scale_y, row=1, col=1)
    fig.update_yaxes(title_text="Inter contact duration distribution", type=scale_y, row=1, col=2)
    fig.update_yaxes(title_text="Number of contacts distribution", type=scale_y, row=2, col=1)
    fig.update_yaxes(title_text="Weight distribution", type=scale_y, row=2, col=2)

    for i, data in enumerate(quantities):
        a = index[i][0]
        b = index[i][1]
        if scale == 'log':
            counts, bins = np.histogram(data, bins=np.logspace(np.log10(min(data - 0.5)), np.log10(max(data + 0.5))),
                                        density=True)

        else:
            counts, bins = np.histogram(data, bins='auto', density=True)
        bins = 0.5 * (bins[:-1] + bins[1:])
        fig.add_trace(go.Histogram(x=bins, y=counts, showlegend=False), row=a, col=b)

    fig.show()


def compare_quantities(quantities_array, label_array, title='Comparison tij data', scale='linear'):
    fig = make_subplots(rows=2, cols=2)
    index = [[1, 1], [1, 2], [2, 1], [2, 2]]
    colors = ['rgb(31, 119, 180)', 'rgb(255, 127, 14)',
              'rgb(44, 160, 44)', 'rgb(214, 39, 40)',
              'rgb(148, 103, 189)', 'rgb(140, 86, 75)',
              'rgb(227, 119, 194)', 'rgb(127, 127, 127)',
              'rgb(188, 189, 34)', 'rgb(23, 190, 207)']
    markers = ['star-triangle-up', 'circle', 'x', 'diamond']

    if scale == 'log':
        scale_x, scale_y = scale, scale

    elif scale == 'linear':
        scale_x, scale_y = scale, scale

    else:
        scale_x, scale_y = 'linear', 'log'

    # Update xaxis properties
    fig.update_xaxes(title_text="Contact duration", type=scale_x, row=1, col=1)
    fig.update_xaxes(title_text="Intercontact duration", type=scale_x, row=1, col=2)
    fig.update_xaxes(title_text="Number of contacts", type=scale_x, row=2, col=1)
    fig.update_xaxes(title_text="weight", type=scale_x, row=2, col=2)

    # Update yaxis properties
    fig.update_yaxes(title_text="Contact duration distribution", type=scale_y, row=1, col=1)
    fig.update_yaxes(title_text="Inter contact duration distribution", type=scale_y, row=1, col=2)
    fig.update_yaxes(title_text="Number of contacts distribution", type=scale_y, row=2, col=1)
    fig.update_yaxes(title_text="Weight distribution", type=scale_y, row=2, col=2)

    for j, data in reversed(list(enumerate(quantities_array))):
        #print("j", j)
        #print("------------")
        #print("data", data)

        data_label = label_array[j]

        for i in range(4):
            a = index[i][0]
            b = index[i][1]
            data = quantities_array[j][i]

            if scale == 'log':
                counts, bins = np.histogram(data, bins=np.logspace(np.log10(np.min(data - 0.5)),
                                                                   np.log10(np.max(data + 0.5))), density=True)

            else:
                counts, bins = np.histogram(data, bins='auto', density=True)

            bins = np.array([(elt + bins[i + 1]) / 2 for i, elt in enumerate(bins[:-1])])
            non_null_index = np.where(counts != 0)[0]
            bins, counts = bins[non_null_index], counts[non_null_index]
            if i==0:
                np.savetxt("contact-duration-ABP-R.5-v.05-noisepi-N1000-100000-bins.txt", bins)
                np.savetxt("contact-duration-ABP-R.5-v.05-noisepi-N1000-100000-counts.txt", counts)
            elif i==1:
                np.savetxt("intercontact-duration-ABP-R.5-v.05-noisepi-N1000-100000-bins.txt", bins)
                np.savetxt("intercontact-duration-ABP-R.5-v.05-noisepi-N1000-100000-counts.txt", counts)
            elif i==2:
                np.savetxt("number-of-contact-ABP-R.5-v.05-noisepi-N1000-100000-bins.txt", bins)
                np.savetxt("number-of-contact-ABP-R.5-v.05-noisepi-N1000-100000-counts.txt", counts)

            else:
                np.savetxt("weight-ABP-R.5-v.05-noisepi-N1000-100000-bins.txt", bins)
                np.savetxt("weight-ABP-R.5-v.05-noisepi-N1000-100000-counts.txt", counts)


            if j == 0:
                if i == 0:
            

                    fig.add_trace(
                        go.Scatter(x=bins, y=counts, mode='lines', marker={'color': colors[j]}, fillcolor=colors[j],
                                   name=data_label), row=a, col=b)
                else:
                   
                    fig.add_trace(
                        go.Scatter(x=bins, y=counts, mode='lines', marker={'color': colors[j]}, fillcolor=colors[j],
                                   name=data_label, showlegend=False), row=a, col=b)

            else:
                if i == 0:
               

                    fig.add_trace(go.Scatter(x=bins, y=counts, marker={'color': colors[j], 'symbol': markers[j - 1]},
                                             name=data_label, mode='markers'), row=a, col=b)
                else:
              

                    fig.add_trace(go.Scatter(x=bins, y=counts, marker={'color': colors[j], 'symbol': markers[j - 1]},
                                             name=data_label, mode='markers', showlegend=False), row=a, col=b)

    fig.show()

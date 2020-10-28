import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
from tqdm import tqdm
import itertools
import matplotlib.mlab as ml
from scipy.interpolate import griddata
from matplotlib import pyplot
from mpl_toolkits.mplot3d import Axes3D
import random
import itertools
import sys

path = "".join(sys.path)
path2 = 'C:\\Users\\Nolan\\SpicePy'
if path2 in path:
    print('in')
else:
    sys.path.append(path2)

import netlist as ntl
from netsolve import net_solve
# read netlist




plt.close('all')

#function to remove the longest edge from the shortest path
def remove_edge_path(path_edges,G):
    edges=[]
    edges_orig=[]
    for things in path_edges:
        info_orig=(G.get_edge_data(things[0],things[1]))
        edges_orig.append(info_orig['weight'])
        if things[0]>3 and things[1]>3:
            info = G.get_edge_data(things[0],things[1])
            edges.append(info['weight'])
        else:
            pass
    edges_np = np.asarray(edges)#convert to np array
    idx = np.argmin(edges_np)#this is the index in the path minus the pins
    idx_edge = edges_orig.index(edges_np[idx])#this is the index in the path
        
    G.remove_edge(path_edges[idx_edge][0],path_edges[idx_edge][1])
    return(G)

#function to build the graph
def build_graph(G,coords,radii_pin,radii_fiber):
    positions={}
    for i in tqdm(range(len(coords))):
        positions[i] = coords[i,:].tolist()
        G.add_node(i,pos=tuple(positions[i]))
        d=(np.sqrt(np.square(coords[i,0]-coords[:,0])+np.square(coords[i,1]-coords[:,1])+np.square(coords[i,2]-coords[:,2])))
        if i<4:
            log = d<radii_pin
        else:
            log = d<radii_fiber
        idxs=np.where(log)[0]
        weights = np.square(d[idxs])
        for j in range(len(idxs)):
            if i!=idxs[j]:
                G.add_edges_from([(i,idxs[j],{'weight':weights[j]})])
    return(G)

#plot the graph
def plot_graph(G,path,path_edges,plot_path=0):

    plt.figure(figsize=(8, 8))
    nx.draw_networkx_edges(G, pos,alpha=0.4)
    
    if plot_path==1:
        nx.draw_networkx_nodes(G,pos,nodelist=path,node_color='r')
        nx.draw_networkx_edges(G,pos,edgelist=path_edges,edge_color='r',width=10)

    plt.xlim(-0.05, 1.05)
    plt.ylim(-0.05, 1.05)
    plt.axis('off')
    plt.show()
    
def connections(pin_idxs,cs):
    poss_conn_all = np.asarray(list(itertools.combinations(pin_idxs,2)))
    poss_conn_undesired=[]
    for i in range(len(poss_conn_all)):
        for j in range(len(cs)-1):
            if poss_conn_all[i,0]==cs[j,0] and poss_conn_all[i,1]==cs[j,1] or poss_conn_all[i,1]==cs[j,0] and poss_conn_all[i,0]==cs[j,1] or poss_conn_all[i,0]==cs[j+1,0] and poss_conn_all[i,1]==cs[j+1,1] or poss_conn_all[i,1]==cs[j+1,0] and poss_conn_all[i,0]==cs[j+1,1]:
                pass
            else:
                poss_conn_undesired.append((poss_conn_all[i,0],poss_conn_all[i,1]))
    return(poss_conn_undesired)

def check_cs(cs,G):
    for j in range(len(cs)):
        if nx.has_path(G,cs[j,0],cs[j,1]):
            stop=0
        else:
            print('broke desired path')
            stop=1
            return(stop)
    return(stop)
    
def pare_network(G_temp,alpha,nodes):
    stop = 0
    G_new=G
    H=G.__class__()
    H.add_nodes_from(G)
    H.add_edges_from(G.edges)
    for i in range((alpha)):
        if nx.has_path(G,nodes[0],nodes[1]):
            path=nx.dijkstra_path(G,nodes[0],nodes[1],weight='weight')
            path_edges = list(zip(path,path[1:]))
            G_new=remove_edge_path(path_edges,G)
            return(G_new,stop)
        else:
            print('no connection')
            stop=1
            return(G_new,stop)
    return(G,stop)
                    

def find_resistance(G,node):
    G1=G.copy()
    weight = 'weight'
    for (u, v, d) in G1.edges(data=True):
       d[weight] = 1/d[weight]                                              
    L=nx.laplacian_matrix(G1,weight='weight').todense()
    M = np.linalg.pinv(L)
    re=[]
    for i in range(len(node)):
        re.append(M[node[i,0],node[i,0]]+M[node[i,1],node[i,1]]-2*M[node[i,1],node[i,0]])
    return(re)
    
def find_resistances(G,node,coords):
    G1=G.copy()
    weight = 'weight'
    for (u, v, d) in G1.edges(data=True):
       d[weight] = 1/d[weight]                                              
    L=nx.laplacian_matrix(G1,weight='weight').todense()
    M = np.linalg.pinv(L)
    log = []
    for i in range(len(coords)):
        log.append(M[node,node]+M[i,i]-2*M[node,i])
    return(log)

def laplacian_resistance_plot(G,coords,pin_indxs,epochs):
    fig,axes = plt.subplots(2,2)
    L=nx.laplacian_matrix(G,nodelist=range(0,len(coords)))
    M = np.linalg.pinv(L.todense())
    xi = np.linspace(0.,1.,100)
    yi = np.linspace(0.,1.,100)
    x=coords[:,0]
    y=coords[:,1]
    fig,axes = plt.subplots(2,2)
    for pins in pin_idxs:
        node=pins
        log=[] 
        for i in range(len(coords)):
            re=M[node,node]+M[i,i]-2*M[i,node]
            log.append(re)
        z=np.asarray(log)
        #zi=scipy.interpolate.griddata(x,y,z,xi,yi,interp='linear')
        zi = griddata((x,y), z, (xi, yi), method='linear')
        axes.flat[pins].contourf(xi,yi,zi)
    
    plt.savefig('re_' + str(epochs) + '.png')


def cool_plot(G_temp):
    G=[]
    G=G_temp
    ncenter=5
    p = dict(nx.single_source_shortest_path_length(G, ncenter))
    
    plt.figure(figsize=(8, 8))
    nx.draw_networkx_edges(G, pos, nodelist=[ncenter], alpha=0.4)
    nx.draw_networkx_nodes(G, pos, nodelist=list(p.keys()),
                           node_size=80,
                           node_color=list(p.values()),
                           cmap=plt.cm.Reds_r)
    
    plt.xlim(-0.05, 1.05)
    plt.ylim(-0.05, 1.05)
    plt.axis('off')
    plt.show() 

def solve_for_v_i(log,bin_num):
    with open("C:\\Users\\Nolan\\Desktop\\vanilla.txt","r+",encoding='utf8') as f:
            data=f.readlines()
        
    combined_data=[]
    for i in range(len(bin_num)):
        ranges=[[0,8],[8,16],[16,24],[24,32]]
        counter=0
        if int(bin_num[i])==0:
            
            for things in range(ranges[i][0],ranges[i][1]):
                combined_data.append(''.join([data[things].strip(),str(i+4),' ', str(counter+1),' ',str(log[things]),'\n']))
                counter+=1
        else:
            for things in range(ranges[i][0],ranges[i][1]):
                combined_data.append(''.join([data[things].strip(),str(1),' ', str(counter+1),' ',str(log[things]),'\n']))
                counter+=1

    
    with open("C:\\Users\\Nolan\\Desktop\\input.net", 'w+',encoding='utf8') as f:
        f.writelines(combined_data) 
        f.writelines(data[32:])

    
    net = ntl.Network('C:/Users/Nolan/Desktop/input.net')  
    try:
        net_solve(net)
    except:
        return(0,0,False)
    
    resistor_num_list = list(range(1,33))
    resistor_list = []
    currents = []
    for items in resistor_num_list:
        resistor_list.append('R' + str(items))
        currents.append(net.get_current('R' + str(items))[0])
        
    resistor_list_dividers = list(range(33,41))
    voltages = []
    for items in resistor_list_dividers:
        voltages.append(net.get_voltage('R' + str(items))[0])
    return(voltages,currents,True)
       

#hyperparameters
radii_pin = 2
radii_fiber = .25
alpha = 2
ratio = .8

#pin locations
pins_x=[0,0,0,0,8,8,8,8,8,8,8,8]
pins_y=[1,3,5,7,1,2,3,4,5,6,7,8]
pins_z=[1,3,5,7,1,2,3,4,5,6,7,8]
pin_idxs=[0,1,2,3,4,5,6,7,8,9,10,11]


#desired connections
#c1 = [0,2]
#c2 = [1,3]
#cs=np.asarray([c1,c2],np.int32)
#in connections
ins = pin_idxs[0:4]
outs = pin_idxs[4:12]

#number of fibers
fibers = 10000
#random fibers
x=np.random.rand(fibers)*max(pins_x)
y=np.random.rand(fibers)*max(pins_x)
z=np.random.rand(fibers)*max(pins_x)
fibers = [x,y,z]
fibers = np.transpose(np.vstack((x,y,z)))


#pins and concat with fibers
pins = np.transpose(np.vstack((np.asarray(pins_x),np.asarray(pins_y),np.asarray(pins_z))))
coords = np.vstack((pins,fibers))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')


ax.scatter(x, y, z, c='r', marker='o')
ax.scatter(pins_x, pins_y, pins_z, c='b', marker='x')

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

#
#instantiate and create graph        
G = nx.Graph()
#build out nodes and edges
G=build_graph(G,coords,radii_pin,radii_fiber)

fig2 = plt.figure(2)
ax = fig2.add_subplot(111)
nx.draw_shell(G)


#get number of edges
N_edges_orig = G.number_of_edges()
N_edges_remaining = round(N_edges_orig*(1-ratio))


#initialize for loop
N_edges_new = G.number_of_edges()
loop = tqdm(total=(N_edges_orig-N_edges_remaining), position = 0)
re_log=[]
large_re_log=[]
#conns=connections(pin_idxs,cs)
poss_conn_all = np.asarray(list(itertools.combinations(pin_idxs,2)))
i=0
stop=0




epochs=1000

error_record=[]

for z in range(epochs):
    #number to pass in
    num_in = (random.randint(1,16))
    truth = num_in*num_in
    bin_num = bin(num_in).split('b')[1]
    bin_num_truth = bin(truth).split('b')[1]
    
    rev=bin_num[::-1]
    num_list = []
    for i in range(0,4):
        try:
            num_list.append(rev[i])
        except:
            num_list.append('0')
    
    num_list_join="".join(num_list)[::-1]
    bin_num=num_list_join
    
    rev=bin_num_truth[::-1]
    num_list = []
    for i in range(0,8):
        try:
            num_list.append(rev[i])
        except:
            num_list.append('0')
    
    num_list_join="".join(num_list)[::-1]
    bin_num_truth=num_list_join
    
    #find the weights of the edges
    #log = find_resistances(G,0,coords)
    weight = 'weight'
    for (u, v, d) in G.edges(data=True):
       d[weight] = 1/d[weight]  
    
    #calculate the resistance
    L=nx.laplacian_matrix(G,weight='weight').todense()
    M = np.linalg.pinv(L)
    
    #calculate all the resistances between pins
    log = []
    connections=[]
    for r in itertools.product(ins, outs):
        connections.append((r[0],r[1]))
    for pairs in connections:
        log.append(M[pairs[0],pairs[0]]+M[pairs[1],pairs[1]]-2*M[pairs[0],pairs[1]])
    
    
    [voltages,currents,test]=solve_for_v_i(log,bin_num)
    if test==True:
        print('worked')
        
        threshold = 4.0
        array=[]
        for things in bin_num_truth:
            array.append(int(things))
        array = np.asarray(array)
        truth_bin_threshold = threshold*array
        voltages_np = np.asarray(voltages)
        error_list = abs(voltages_np-truth_bin_threshold)
        error = ((voltages_np-truth_bin_threshold)**2).mean(axis=0)
        v_idx_to_prune = voltages.index(max(voltages))
        c_idx_to_prune = currents.index(max(currents))
        out_pin_to_prune=outs[v_idx_to_prune]
        
        idx_list= []
        for i in range(len(connections)):
            things = connections[i]
            if things[1] == out_pin_to_prune:
                 idx_list.append([things[0],things[1],i]) 
        
        currents_to_test=[]
        for item in idx_list:
            currents_to_test.append(currents[item[2]]) 
        max_current_idx = currents_to_test.index(max(currents_to_test))
        connections_to_prune = idx_list[max_current_idx][0:2]
        
        
        [test,stop]=pare_network(G,alpha,connections_to_prune)
        if stop==1:
            break
        
        error_record.append(error)
        
fig = plt.figure(3)
ax = fig.add_subplot(111)
ax.plot(error_record)
plt.xlabel('Epochs')
plt.ylabel('error')
plt.title('Training error')
#loop while threshold has not been reached (number of edges removed)
#while N_edges_new>N_edges_remaining:
#    N_edges_old = N_edges_new
#    #pare the network
#    [G_new,stop]=pare_network(conns,G,alpha,cs)
#    
#    #stop the training if the stop criteria is met
#    if stop==1:
#        print('stopping')
#        #cool_plot(G_new)
#        break
#    
#    #find edge-removal progress
#    N_edges_new = G.number_of_edges()
#    delta_edges = N_edges_old-N_edges_new
#    loop.update(delta_edges)
#    
#    #get the shortest connections for all of the pins
#    temp=[]
#    for j in range(len(poss_conn_all)):
#        if nx.has_path(G,poss_conn_all[j,0],poss_conn_all[j,1]):
#            temp.append(nx.dijkstra_path_length(G,poss_conn_all[j,0],poss_conn_all[j,1]))
#        else:
#            temp.append(float("inf"))
#    if i==0:
#        large_re_log=temp
#    else:
#        large_re_log=np.vstack((large_re_log,temp))
#    if i%10==0:
#        if i==0:
#            re_real_log=find_resistance(G,poss_conn_all)
#        else:
#            re_real_log=np.vstack((re_real_log,find_resistance(G,poss_conn_all)))
#    i+=1
#
#    #if i%50==0:
#    #    laplacian_resistance_plot(G,coords,pin_idxs,i)
##    if i%20==0:
##        re=find_resistance(G,cs[0,:])
##        re_log.append(re)           
#        
#plt.figure(0)
#plt.plot(large_re_log)
#plt.legend(('1','2','3','4','5','6'))
#plt.figure(1)
#plt.plot(re_real_log)
#plt.legend(('1','2','3','4','5','6'))
#
#
#
#
##laplacian_resistance_plot(G,coords,pin_idxs,i)
#
#
#log = find_resistances(G,0,coords)
#
#fig = pyplot.figure()
#ax = Axes3D(fig)
#norm = plt.Normalize()
#log_test=np.square(np.square(np.asarray(log)))
#colors = plt.cm.jet(norm(log_test))
#
#ax.scatter(coords[:,0], coords[:,1], coords[:,2], c = colors)
#pyplot.show()
#
#path=nx.dijkstra_path(G,1,2,weight='weight')
#path_coordsx=coords[path,0]
#path_coordsy=coords[path,1]
#path_coordsz=coords[path,2]
#
#
#ax.scatter(path_coordsx, path_coordsy, path_coordsz, c = 'r')
#pyplot.show()
#
#fig2 = plt.figure(10)
#ax = Axes3D(fig2)
#ax.scatter(path_coordsx, path_coordsy, path_coordsz, c = 'r')
#pyplot.show()
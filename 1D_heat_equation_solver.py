import numpy as np
import matplotlib.pyplot as plt

k=x=dx=t=dt=grid_rows=grid_columns=length_nodes=time_nodes=alpha=T1=T2=t_initial=0
grid=left_matrix=[]

def get_input(method):
    '''This function gets inputs required from the user'''

    global k,x,dx,t,dt
    if(method.upper()=='A' or method.upper()=='C'):
        k=float(input('Enter diffusivity(k)(in cm^2/s): '))
        x=float(input('Enter length of rod(in cm): '))
        dx=float(input('Enter length step size: '))
        t=float(input('Enter simulation time(in s): '))
        dt=float(input('Enter time step size: '))
        if ((dx<=0 or dt<=0) or (dx>x or dt>t)):
            raise ValueError

    elif(method.upper()=='B' or method.upper()=='D'):
        k=float(input('Enter diffusion coefficient(in cm^2/s): '))
        x=float(input('Enter length of rod(in cm): '))
        t=float(input('Enter simulation time(in s): '))



def initialize_grid():
        ''' This function initializes the length-time grid with zeros'''

        global grid_rows,grid_columns,time_nodes,length_nodes,grid

        grid_rows=time_nodes+1
        grid_columns=length_nodes+2
        grid=np.zeros((grid_rows,grid_columns))



def boundary_conditions():
        '''This function gets initial and boundary conditions from the user and updates it in the grid'''

        global T1,T2,grid,grid_rows,length_nodes,time_nodes,t_initial
        print('Enter initial and boundary conditions(in celcius): ')
        T1=int(input('Enter surface 1 temperature: '))
        T2=int(input('Enter surface 2 temperature: '))

        for i in range(grid_rows):
            grid[i,0]=T1
            grid[i,length_nodes+1]=T2
        t_initial=int(input('Enter initial rod temperature: '))
        for i in range(1,length_nodes+1):
            grid[time_nodes,i]=t_initial
        #print(grid)



def calculate_finite_difference_method():
    '''This function takes in the grid with boundary conditions applied and calculates temperature
    at other nodes using finite difference method and updates it to the grid '''

    global time_nodes,length_nodes,alpha,grid
            #first loop increments row of the grid(i.e)Time
            #second loop increments column of the grid(i.e)length

    for i in range(time_nodes,0,-1):
            for j in range(1,length_nodes+1):
                U=grid[i,j]                             #temperature at current node
                Uw=grid[i,j-1]                          #temperature at node to the west
                Ue=grid[i,j+1]                          #temperature at node to the east

                Un=U+(alpha*(Ue-2*U+Uw))            #calculates temperature at current node at time t+1
                grid[i-1,j]=round(Un,4)                 

    #print('grid',grid)



def generate_tri_diagonal_matrix():
    '''This function generates a tri diagonal matrix using the value of alpha '''

    global alpha,length_nodes,left_matrix
    diag=2*(1+alpha)
    below_diag=-alpha
    above_diag=-alpha

    left_matrix=np.zeros((length_nodes,length_nodes))
    for i in range(length_nodes-1):
                left_matrix[i][i]=diag
                left_matrix[i][i+1]=above_diag
                left_matrix[i+1][i]=below_diag
    left_matrix[length_nodes-1][length_nodes-1]=diag

def calculate_crank_nicolson_method():
    '''This function calculates the temperature of all nodes at time:t+1 using left and coefficient matrix and updates it to the grid'''
    
    global length_nodes,time_nodes,grid,left_matrix

    temp=np.zeros((length_nodes,1))                         #temperature matrix(unknown)
    coeff_matrix=np.zeros((length_nodes,1))                 

    #based on the position of node the coefficient matrix is calculated

    for i in range(time_nodes,0,-1):
        for j in range(1,length_nodes+1,1):
            #left node
            if j==1:
                U=grid[i,j]                                     
                Uw=grid[i,j-1]                                  
                Ue=grid[i,j+1]                                  
                Unw=grid[i-1,j-1]                               #temperature at node to the north-west
                Une=grid[i-1,j+1]                               #temperature at node to the north-east
                coeff_matrix[j-1,0]=(alpha*(Uw+Unw+Ue))+(2*(1-alpha)*U)    
            
        #right node
            elif j==(length_nodes):
                U=grid[i,j]
                Uw=grid[i,j-1]
                Ue=grid[i,j+1]
                Une=grid[i-1,j+1]
                coeff_matrix[j-1,0]=(alpha*(Uw+Ue+Une))+(2*(1-alpha)*U)
            
            
            #center node
            else:
                U=grid[i,j]
                Uw=grid[i,j-1]
                Ue=grid[i,j+1]
                coeff_matrix[j-1,0]=(alpha*(Uw+Ue))+(2*(1-alpha)*U)
        
    #Similar to AX=B here we have left_matrix*temp=coeff_matrix
        temp= np.linalg.inv(left_matrix).dot(coeff_matrix)
        for k in range(1,length_nodes+1):
                grid[i-1,k]=round(temp[k-1,0],2)                



def display_output():
    '''This function displays output in two forms:
    1. Temperature at every node on the rod at the end of simulation time
    2.Graph showing how the temperature has changed over simulation time w.r.t time and length'''

    global grid_columns,grid,time_nodes,grid_rows,T1,T2,length_nodes,dx,dt

    print('Temperature of rod at the end of simulation time: ')
    for i in range(grid_columns):
        print(grid[0,i],end=" ")

    #Graph shows 5 plots with time difference=simulation time/5

    if time_nodes>5:
        X=[]                        
        Y=[]
        d=round(grid_rows/5)                 #divides total time into 5                    
        time=0
        for i in range(time_nodes,-1,-d):
            X.append(0)                                         
            Y.append(T1)
        
            for j in range(1,length_nodes+1): 
                X.append(j*dx)
                Y.append(grid[i,j])

            X.append((length_nodes+1)*dx) 
            Y.append(T2)
        
            s="at t ="+str(time)+"s"
            plt.plot(X,Y,label=s) 
            plt.xlabel("Bar length")
            plt.ylabel("Temperature")
            time=round((time+(dt*d)),2)
            X.clear()
            Y.clear()
        plt.legend()
        plt.show()


def main():
    global k,x,dx,t,dt,grid_rows,grid_columns,length_nodes,time_nodes,alpha,T1,T2,t_initial,grid,left_matrix
    
    method=input("Select a method to solve 1D heat equation\n a.finite difference method \n b.bender schmidt method(λ=0.5) \n c.crank nicolson method \n d.crank nicolson method(λ=1) \n")

    if(method.upper()=='A'):
        get_input('A')
        length_nodes=int(x//dx)-1
        time_nodes=int(t//dt)
        alpha= k*dt/(dx**2)
        if alpha<=0.5:
            initialize_grid()
            boundary_conditions()
            calculate_finite_difference_method()
            display_output()
        else:
            print('alpha value exceeds 0.5 ,can not use forward difference difference method')

    if(method.upper()=='B'):
        get_input('B')
        dx=1
        dt=0.5*k
        length_nodes=int(x//dx)-1
        time_nodes=int(t//dt)
        alpha= 0.5
        initialize_grid()
        boundary_conditions()
        calculate_finite_difference_method()
        display_output()

    if(method.upper()=='C'):
        get_input('C')
        length_nodes=int(x//dx)-1
        time_nodes=int(t//dt)
        alpha= k*dt/(dx**2)
        initialize_grid()
        boundary_conditions()
        generate_tri_diagonal_matrix()
        calculate_crank_nicolson_method()
        display_output()


    if(method.upper()=='D'):
        get_input('D')
        dx=1
        dt=1/k
        length_nodes=int(x//dx)-1
        time_nodes=int(t//dt)
        alpha=1

        initialize_grid()
        boundary_conditions()
        generate_tri_diagonal_matrix()
        calculate_crank_nicolson_method()
        display_output()

if __name__=="__main__":
    main()
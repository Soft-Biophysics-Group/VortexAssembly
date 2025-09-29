#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  6 19:10:42 2023

@author: lara
"""

import json
from SideClasses import SideFunctions as SF
import numpy as np 

import matplotlib 
from matplotlib import cm
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

from matplotlib.patches import FancyArrowPatch




class ReadLatticeParticle():
    
    
    def __init__(self, name , directory):
        
        file = open(directory + name +'.json', 'r')
        
        myDict = json.load(file)
        
        self.name = name
        self.directory=directory
        self.generators = myDict['Generators'] 
        self.edges = np.array(myDict['Edges'])
        self.n_edges = myDict['N_Edges']
        
        self.generators = np.array(self.generators)
        
        self.neighbors =np.array(myDict['Neighbors'])
        self.n_neighbors = myDict['N_Neighbors']
        
        
        self.neighbors_cartesian = np.zeros((self.n_neighbors,3))
        for i in range (self.n_neighbors):
            self.neighbors_cartesian[i] = self.lattice_to_cartesian(self.neighbors[i])
            norm = np.linalg.norm(self.neighbors_cartesian[i])
            self.neighbors_cartesian[i] = self.neighbors_cartesian[i]/norm
        
        self.n_orientations = myDict['N_Orientations']
        self.orientations = np.array(myDict['Orientations'])
        
        
        
        """Usually the only unique edge is edge=1, but in the case of the hexagon 
        we also need edge=4 to get all the structures"""
        self.face_type = myDict['Faces_type']
        self.unique_faces = np.unique(self.face_type)
        self.unique_edges = np.array([np.where(self.face_type==unique)[0][0] for unique in self.unique_faces])+1
        
                                
        
        """Total number of structures"""
        self.N_particles_type_max = 0
        self.all_N_structures = {}
        self.detail_all_N_structures = {}
        maxfound = False 
        i = 1
        while not maxfound:
            try :
                this_n_structures = myDict['N_Structures_'+str(i)]
                detail_this_n_structures = myDict['detail_N_Structures_'+str(i)]
            except :
                maxfound=True
            
            if not maxfound : 
                self.all_N_structures[i] = this_n_structures
                self.detail_all_N_structures[i] = detail_this_n_structures
                self.N_particles_type_max +=1
            i+=1
            
        self.read_structure_file()
        self.build_state_dict()
        self.compute_opposite_orientations()
        
        if name=='Triangular':
            self.map_old_to_new_empty_full()
                
    def read_structure_file(self):
        self.config_to_structure_dicts = {} # (o1, o2, edge) -> s
        self.structure_to_config_dicts = {} # s -> (o1, o2, edge=1), o1<o2
        self.structure_degeneracy_dicts = {} #how many configuratoin have the same s
        
        for ntype in range(1, self.N_particles_type_max):
            self.config_to_structure_dicts[ntype] = {}
            self.structure_to_config_dicts[ntype] = {}
            self.structure_degeneracy_dicts[ntype] = {}
        
            all_structures = SF.matrixToArray(self.directory+self.name+'_structures_'+str(ntype)+'particles.txt')
            for a_structure in all_structures :
                o1, o2, edge, s = a_structure
                
                self.config_to_structure_dicts[ntype][(o1, o2, edge)]=s
      
                if s in self.structure_degeneracy_dicts[ntype].keys():
                    self.structure_degeneracy_dicts[ntype][s]+=1  
                else :
                    self.structure_degeneracy_dicts[ntype][s] = 0
                    

                # if self.name=='Triangular' and s==71 and ntype==2:
                #     print('o1,o2,edge', o1, o2, edge, 's',s)
                        
                if (edge in self.unique_edges  or edge==0) :
                    
                    if not s in self.structure_to_config_dicts[ntype].keys():  
                        self.structure_to_config_dicts[ntype][s]=(o1, o2, edge)
                        
    def get_structure(self, n_particles_type,  o1,o2, edge, twoParticlesStructures=True):
        """Returns the indices as saved for python and not for c++, without the initial one particle structures"""
        N1 = self.detail_all_N_structures[n_particles_type][0]
        return(self.config_to_structure_dicts[n_particles_type][(o1, o2, edge)]-N1)
                    
    def get_particle_state(self, particle, orientation):
        if orientation==0:
            return(0)
        else :
            return(particle*self.n_orientations+orientation)
    
    def get_particle_id_orientation(self, particle_state):
        if particle_state ==0 :
            particle_id=0
            orientation=0
        else :
            particle_id = (particle_state-1)//self.n_orientations
            orientation = (particle_state-1)%self.n_orientations+1
        return(particle_id, orientation)
                    
    def build_state_dict(self):
        self.particle_orientation_to_state = {}
        self.state_to_particle_orientation = {}
        for particle in range (self.N_particles_type_max):
            for orientation in range (self.n_orientations+1):
                state = self.get_particle_state(particle, orientation)
                self.particle_orientation_to_state[(particle, orientation)] = state
                self.state_to_particle_orientation[state] = (particle, orientation)
        
        
    def lattice_to_cartesian(self, coordinate):
        """Transforms lattice coordinate to cartesian coordinate"""
        d =  len(self.generators)
        x,y,z = 0,0,0
        for i in range(d):
            dx, dy, dz = coordinate[i] * self.generators[i]
            x+=dx
            y+=dy
            z+=dz
        return(x,y,z)
    
    

    def lattice_resquare(self, coordinate, Lx, Ly):
        """In lattice coordinates"""
        x,y,z = coordinate
        y = y%Ly
        D=  y//2
        
        if x>Lx - D -1 :
            x = x-Lx
        if x<-D:
            x= x + Lx 
        return(np.array([x,y,z]))
    
    def find_opposite_edge_id(self, edge_id):
        """If edge_id is the id of the vector pointing towards the particle
        what is the id of the same vector pointing outwards the particle"""

        if edge_id ==0 :
            return(0)
        else :
            edge = self.neighbors[edge_id-1]
            opposite_edge = -edge
            searchZeros = np.abs(self.neighbors -opposite_edge)
            opposite_edge_id = np.where((searchZeros==0).all(axis=1))[0][0]
            opposite_edge_id +=1
            return(opposite_edge_id)


    def find_face_id(self, orientation, edge_id):
            """If particle is in given orientation, 
            what face is facing the edge number edge_id
            If the orientation is 0, then the face is -1"""
            if orientation==0 :
               face_id = 0
            else :
                face_id = np.where(self.orientations[orientation-1] == edge_id-1)[0][0]+1
            return(face_id)
    
    
    def get_faces_in_contact(self, structure, ntype):
        """Returns particle1, face1, particle2, face2"""        
        state1, state2, edge = self.structure_to_config_dicts[ntype][structure]  
        if edge!=0 :
            if state1!=0 :
                p1, o1 = self.state_to_particle_orientation[state1]
                face1 = self.find_face_id(o1, edge) 
            else :
                p1 = 0
                face1 = 0
            if state2 !=0:
                p2, o2 = self.state_to_particle_orientation[state2]               
                face2 = self.find_face_id(o2, self.find_opposite_edge_id(edge))
            else :
                p2 = 0
                face2 = 0
        else :
            p1, face1, p2, face2 = 0,0,0,0
        
        return(p1, face1, p2, face2)
    
    
    def show_structure(self, s=None, 
                       state1=None, state2=None, edge=None,
                       p1=None, p2=None,o1 = None, o2=None, 
                       ntype=1, axes=None, representation='arrowed', dimension=2):
        """Structure is either defined by s, by (state1, state2, edge), or by (p1, o1, p2, o2, edge)"""
        
        if o1 is not None :
            state1 = self.get_particle_state(p1, o1)
            state2 = self.get_particle_state(p2, o2)
            
            # print(state1, state2)
        
        if state1 is not None :
            s = self.config_to_structure_dicts[ntype][(state1, state2, edge)]
            
        
        state1, state2, edge = self.structure_to_config_dicts[ntype][s]
    
        p1, o1 = self.state_to_particle_orientation[state1]
        p2, o2 = self.state_to_particle_orientation[state2]
            
        
        
        print('p1,p2,o1,o2,edge', p1,p2,o1,o2,edge)
        size=0.95
        Particle = ParticleRepresentation(self)
        x1, y1, z1 = 0, 0, 0
        dx, dy, dz = self.neighbors_cartesian[edge-1]
        x2, y2, z2 = x1+dx, y1+dy, z1+dz
        # print(x1, x2, y1, y2)
        if axes is None :
            f,axes=plt.subplots(figsize=(3,3))
        print(axes, representation)
        Particle.plot_particle(x1, y1,z1,  size, p1, o1, dimension, axes=axes, representation=representation)
        Particle.plot_particle(x2, y2, z2, size, p2, o2, dimension, axes=axes, representation=representation)
        
        xmin, xmax =min(x1,x2), max(x1,x2)
        ymin, ymax =min(y1,y2), max(y1,y2)
        axes.set_xlim(xmin-size, xmax+size)
        axes.set_ylim(ymin-size, ymax+size)
        axes.set_title('s='+str(s))
        axes.set_aspect('equal', 'box')  
        axes.tick_params(axis='both',  bottom=False, left=False, labelbottom=False, labelleft=False)
        return(axes,s )
    
    
    
    def compute_opposite_orientations(self):
        """Creates a dict that for orientation o1 gives the reverse orientation """
        orientations = np.arange(self.n_orientations+1)
        opposites = - np.ones(self.n_orientations+1)
        for o1 in orientations :
            for o1op in range(o1+1) :
                s1 = self.config_to_structure_dicts[1][o1, 0, 1]
                s2 = self.config_to_structure_dicts[1][0, o1op, 1]
                if s1==s2 :
                    opposites[o1] = o1op
                    opposites[o1op] = o1
        
        self.opposite_orientations = opposites
    
    
    
    def get_matrix_structures(self, p1=0, p2=0, symmetric=True, n_type=1, emptyFullAlso=False):
        n = self.n_orientations + emptyFullAlso
        
        o1s = np.arange(1, n+1)
        if symmetric :
            o2s = self.opposite_orientations.copy()[1:]
        else :
            o2s = o1s
            
        if emptyFullAlso:
            o1s = np.concatenate((np.zeros(1), o1s))
            o2s = np.concatenate((np.zeros(1), o2s))
        
        structures_matrix = np.zeros((n,n), dtype=int)
        for i1 in range (n):
            for i2 in range (n):
                o1, o2 = o1s[i1], o2s[i2]
                state1, state2 = self.particle_orientation_to_state[p1, o1], self.particle_orientation_to_state[p2, o2]
                structures_matrix[i1, i2] = self.config_to_structure_dicts[n_type][state1, state2, 1]
        return(structures_matrix)
    
    
    def map_old_to_new_empty_full(self):
        matrice_structure = self.get_matrix_structures()
        
        faces_order_old = [3,2,1,6,5,4]
        # faces_order_old = [2,3,4,5,6,1]
        face_to_empty_structure_old = {}
        face_to_empty_structure_old[1]=6
        face_to_empty_structure_old[6]=3
        face_to_empty_structure_old[5]=2
        face_to_empty_structure_old[4]=1
        face_to_empty_structure_old[3]=4
        face_to_empty_structure_old[2]=5
        
        face_to_empty_structure_new = {}
        for o in range (1,7):
            s = self.config_to_structure_dicts[1][o,0,1]
            p,f,p,zero = self.get_faces_in_contact(s,1)
            face_to_empty_structure_new[f]= s
        
        self.old_to_new = {}
        self.new_to_old = {}
        
        for i in range(6):
            """which faces are in contact in the ith diagonal item of the ref matrix"""
            s_diag = matrice_structure[i][i]
            p1, face, p2, face = self.get_faces_in_contact(s_diag, 1)
            """To which empty-full strucutre this face correspond"""
            s_empty = face_to_empty_structure_new[face] -7
            
            """in the old matrix, the ith diagonal index correspond to old_face"""
            old_face = faces_order_old[i]
            """To which empty-full strucutre  old_face correspond"""
            old_s_empty = face_to_empty_structure_old[old_face]
            
            self.old_to_new[old_s_empty] = s_empty
            self.new_to_old[s_empty] = old_s_empty
    

       
            



    def full_structures_each_face(self, n_structures_type=1):
        s_faces = np.zeros((self.n_neighbors, self.n_neighbors), dtype=int)
        for s in range(self.all_N_structures[n_structures_type]):
            p1, f1, p2, f2 = self.get_faces_in_contact(s, n_structures_type)
            if f1*f2 !=0 :
                s_faces[f1-1][f2-1] = s
                s_faces[f2-1][f1-1] = s
        return(s_faces)
    
    
    
    def rotate_local_structure(self, local_structure, n):
        
        o, s1, s2, s3, s4, s5, s6 = local_structure
        structures = [s1, s2, s3, s4, s5, s6]
        
        o_new = (o-n-1)%6+1
        structures_new = structures[n:]+structures[:n]
        
        return([o_new]+structures_new)
    
    def local_structure_to_ref(self, local_structure):
        
        o, s1, s2, s3, s4, s5, s6 = local_structure
        n = (o-1 -1)%6+1
        return(self.rotate_local_structure(local_structure, n))
        
            
        
                    
class ParticleRepresentation():
    """
    This class is used to plot particles in 2D or 3D knowing the lattice they belong to
    You can have several choice for the representation, which is then an argument of the plotting function
    To plot a particle, you need to know its position, orientation, and type of particle
    """
    
    
    
    
    def __init__(self, lattice ):
        self.lattice = lattice
        self.n_faces = self.lattice.n_neighbors
        self.orientations = self.lattice.orientations
        
        self.contact_centers_polar = [2*np.pi*(k/self.n_faces) for k in range(self.n_faces)]    
        self.angular_increment =  self.contact_centers_polar[1]
        # the theta coordinate for the centers of the face of the polyhedra
        self.vertices_polar = np.array(self.contact_centers_polar) +self.angular_increment/2
        # the theta coordinate for the vertices of the face of the polyhedra
        
        self.cos_contacts = np.cos(self.contact_centers_polar)       
        self.sin_contacts = np.sin(self.contact_centers_polar)      
        self.cos_vertices = np.cos(self.vertices_polar)
        self.sin_vertices = np.sin(self.vertices_polar)
        
        
        
        
        
        
        """Cube vertices"""
        #front face 
        # 3   4
        # 1   2
        # back face :
        # 7   8
        # 5   6
        v1 = np.array([0,0,0])
        v2 = np.array([1,0,0])
        v3 = np.array([0,0,1])
        v4 = np.array([1,0,1])
        v5 = np.array([0,1,0])
        v6 = np.array([1,1,0])
        v7 = np.array([0,1,1])
        v8 = np.array([1,1,1])
        self.cubeVerticies = [v1, v2, v3, v4, v5, v6, v7, v8]
        
        
        # colormaps = [cm.viridis, cm.seismic]
        colormaps = [cm.viridis, cm.spring, cm.summer, cm.winter, cm.autumn,cm.cool]
        imax = [0.9, 0.8, 0.8, 0.8, 0.5, 0.8]
        self.myColors1 = [colormaps[i](np.linspace(0,imax[i], self.n_faces)) for i in range(len(colormaps))]
        
        self.myColors2 = ['gray', 'silver', 'royalblue', 'mediumorchid']
        self.myColors2 = ['snow', 'gold', 'mediumaquamarine', 'royalblue', 'mediumorchid']
        
        """FFL, LLK"""
        self.myColors3 = [[ 'dodgerblue','dodgerblue','gold','gold','gold','gold'], 
                          ['darkblue', 'darkblue','dodgerblue', 'dodgerblue','dodgerblue', 'dodgerblue']]
        self.myColors3 = [['darkblue','darkblue','dodgerblue','dodgerblue','gold','gold'], 
                          ['darkblue', 'darkblue','dodgerblue', 'dodgerblue','dodgerblue', 'dodgerblue']]
        
        
        
        self.myColorPatchy = ['royalblue', 'gold', 'salmon', 'forestgreen', 'rebeccapurple', 'darkorange']
        
        self.add_angle = self.n_faces==4
        
    def plot_2Dtriangular_patch(self,x0, y0, size,edge_id, color_id, axes):
        
        ratio_size = 0.5/self.cos_vertices[0]
        
        phi = self.contact_centers_polar[edge_id] -self.angular_increment/2
        
        path = [[x0,y0]]
        i1, i2 = edge_id, (edge_id+1)%self.n_faces
        path += [[x0+size*ratio_size*self.cos_vertices[i1], y0+size*ratio_size*self.sin_vertices[i1]]]
        path += [[x0+size*ratio_size*self.cos_vertices[i2], y0+size*ratio_size*self.sin_vertices[i2]]]
        path=np.array(path)
        
        h = mpatches.Polygon(path, color=self.myColors1[color_id[0]][color_id[1]], linewidth=0)#,ec='black')
        axes.add_artist(h)
        
    def plot_2Dtriangular_patch_triangular(self, x0, y0, size, edge_id, color_id, axes ):
        
        ratio_size = 0.5/self.cos_vertices[0] *1.5
        
        phi = self.contact_centers_polar[edge_id] -self.angular_increment/2
        
        path = [[x0,y0]]
        i1, i2 = (edge_id)%self.n_faces, (edge_id+2)%self.n_faces
        path += [[x0+size*ratio_size*self.cos_contacts[i1], y0+size*ratio_size*self.sin_contacts[i1]]]
        path += [[x0+size*ratio_size*self.cos_contacts[i2], y0+size*ratio_size*self.sin_contacts[i2]]]
        path=np.array(path)
        
        h = mpatches.Polygon(path, color=self.myColors3[color_id[0]][color_id[1]], linewidth=0)#,ec='black')
        axes.add_artist(h)
        
        
        
        
    def plot_2Dtriangle(self, x0, y0, size, particle, orientation, axes, order=1):
        if orientation > 0 :
            if orientation%2 ==0:
                i_edges = [0,2,4]
                
            else :
                i_edges = [1,3,5]
                
            for i_edge in i_edges :
                color_face = (self.orientations[1-orientation][i_edge]* order)%self.n_faces
                color_id = (particle, color_face)
                self.plot_2Dtriangular_patch_triangular(x0, y0, size,i_edge, color_id, axes)
        
        
        
        
        
    def plot_2Dcolored_particle(self, x0, y0, size, particle, orientation, axes, order=1):
        if orientation > 0 :
            for i_edge in range (self.n_faces):
                color_face = (self.orientations[1-orientation][i_edge]* order)%self.n_faces
                color_id = (particle, color_face)
                self.plot_2Dtriangular_patch(x0, y0, size,i_edge, color_id, axes)
                
                
    def plot_2Dpatchy_particle(self, x0, y0, size, particle, orientation, axes, patchID=[0,1,2,3,4,5], dotsize=100):
        
        
        if orientation >0 :
            
            radius = 0.5/self.cos_vertices[0] * size
            ratio_size = 3#0.5/self.cos_vertices[0] *1.5
            phi = self.angular_increment/2 *self.add_angle
            h = mpatches.RegularPolygon((x0,y0), self.n_faces, radius=radius, orientation=phi, facecolor="silver", edgecolor='black')
            axes.add_artist(h)
            
            for i_edge in range (self.n_faces):
                color_id = patchID[self.orientations[1-orientation][i_edge]]
                x = x0+radius*self.cos_vertices[i_edge]
                y= y0+radius*self.sin_vertices[i_edge]
                axes.scatter(x,y,s=dotsize, color=self.myColorPatchy[color_id], zorder=1e4)
            
           
        
    
    def plot_2Darrowed_particle(self, x0,y0,size, particle, orientation, axes, edgecolor='none'):
        if orientation > 0 :
            radius = 0.5/self.cos_vertices[0] * size
            phi = self.angular_increment/2 *self.add_angle
            h = mpatches.RegularPolygon((x0,y0), self.n_faces, radius=radius, 
                                        orientation=phi, facecolor=self.myColors2[particle],
                                        edgecolor=edgecolor)
            
            arrow_radius = 1*radius/2
            
            
            face_id = self.orientations[orientation-1][0] 
            xh,yh = self.cos_contacts[face_id],self.sin_contacts[face_id]
            xt, yt = -xh, -yh
            xa0, ya0 = x0+xt*arrow_radius , y0+yt*arrow_radius
            dxa, dya = (xh-xt)*arrow_radius, (yh-yt)*arrow_radius
            
            head_length = arrow_radius/2
            head_width= arrow_radius/3
            width=size/40
            
            arrow = mpatches.FancyArrow(xa0, ya0, dxa, dya, width=width, color='black', length_includes_head=True, head_width=head_width,head_length=head_length)
                
            axes.add_artist(h)
            axes.add_artist(arrow)
        
        if orientation==0:
            radius = 0.5/self.cos_vertices[0] * size
            phi = self.angular_increment/2 *self.add_angle
            h = mpatches.RegularPolygon((x0,y0), self.n_faces, radius=radius, 
                                        orientation=phi, facecolor='white',
                                        edgecolor=edgecolor)
            axes.add_artist(h)
            
            
    def plot_2Dparticle_orientation(self, x0, y0, size, particle, orientation, axes, order):
        if orientation>0:
            radius = 0.5/self.cos_vertices[0] * size
            phi = self.angular_increment/2 *self.add_angle
            h = mpatches.RegularPolygon((x0,y0), self.n_faces, radius=radius, 
                                        orientation=phi, facecolor=self.myColors1[particle][orientation-1])
                                        # ,edgecolor=edgecolor)
            axes.add_artist(h)
            
        
    def plot_2Dparticle(self, x0, y0, size, particle, orientation, axes, representation, order=1):
        if representation=="colored":
            self.plot_2Dcolored_particle(x0, y0, size, particle, orientation, axes, order=order)
        if representation=="arrowed":
            self.plot_2Darrowed_particle( x0,y0,size, particle, orientation, axes)
        if representation=="triangle":
            self.plot_2Dtriangle( x0,y0,size, particle, orientation, axes)
        if representation=="orientation":
            self.plot_2Dparticle_orientation(x0, y0, size, particle, orientation, axes, order=order)
            


    def plot_3Darrow(self, x0, y0, z0, size, particle, orientation, axes):
        if orientation >0 :  
            
            face1, face2, face3 = 1,2,3
            edges = []
            for face in [face1, face2, face3]:
                edges.append(self.orientations[orientation-1][face-1] + 1)
            edge1, edge2, edge3 = edges
            
            color = self.myColors2[particle]
            
            lw_0 = 1
            lw_ratio = 1
            for edge in edges[:2] :
                dx, dy, dz = self.lattice.neighbors_cartesian[edge-1]
                
                arrow_radius = size/lw_ratio/2
                x1,x2 = x0-dx*arrow_radius, x0+dx*arrow_radius
                y1,y2 = y0-dy*arrow_radius, y0+dy*arrow_radius
                z1,z2 = z0-dz*arrow_radius, z0+dz*arrow_radius
            
                lw = lw_0/lw_ratio
                lw_ratio+=1
                mutation_scale = 3*size
                # head_width = arrow_radius/3
                # head_length = arrow_radius/2
                arrow = myArrow3D([x1, x2], [y1,y2], [z1,z2], mutation_scale=mutation_scale, lw=lw,
                                  arrowstyle="simple", color=color)
                                  # , length_includes_head=True, 
                                  # head_width=head_width,head_length=head_length)
                axes.add_artist(arrow)


    def plot_Cube(self, x0, y0, z0, size, particle,  orientation , axes, alpha=0.95):
        """alpha is the transparency """
        if orientation > 0:
            
            """Get the colors of all the faces"""
            myEdges = self.lattice.orientations[orientation-1]
            faces_colors = []
            for edge in range(6):
                # print(myEdges, edge, myEdges-edge)
                face = np.where(np.array(myEdges)==edge)[0][0]
                faces_colors.append(self.myColors1[particle][face])
            
            """Get the coordinate of the vertices"""
            radius = size/2
            center= np.array([x0,y0,z0])
            vertices = self.cubeVerticies
            new_vertices = []
            for v in vertices :
                new_vertices.append( center-radius + 2*radius*v)
            v1, v2, v3, v4, v5, v6, v7, v8 = new_vertices
            
            face_1 = np.array([[v2,v6],[v4,v8]])
            face_2 = np.array([[v6,v5],[v8,v7]])
            face_3 = np.array([[v4,v8],[v3,v7]])
            
            face_4 = np.array([[v1,v5],[v3,v7]])
            face_5 = np.array([[v1,v2],[v3,v4]])
            face_6 = np.array([[v1,v2],[v5,v6]])
            
            faces = [face_1, face_2, face_3, face_4, face_5, face_6]
            
            """Plot each face"""
            for i in range(len(faces)):
                face = faces[i]
                X,Y,Z = face[:,:,0], face[:,:,1], face[:,:,2]
                color = faces_colors[i]
                color[3] = alpha
                axes.plot_surface(X, Y, Z ,color=color, edgecolor='none')
            
        
    def plot_3Dhexagon(self, x0, y0, z0, size, particle, orientation, axes, representation=None):
        """Plots a 3D hexagon in a 2d lattice, with deferent colors depending on if it's up or donw"""
        
        orientation_corresp = {}
        orientation_corresp[10]=1
        orientation_corresp[9]=2
        orientation_corresp[11]=3
        orientation_corresp[8]=4
        orientation_corresp[12]=5
        orientation_corresp[7]=6 
        Triangular = ReadLatticeParticle('Triangular', '/Users/lara/Documents/SelfAssembly2/Lattice/')
        newRepresentation = ParticleRepresentation(Triangular)
        if orientation in [1,2,3,4,5,6]:            
            newRepresentation.plot_2Dparticle(x0, y0, size, particle, orientation, axes, representation='arrowed')
        
        else :
            newRepresentation.plot_2Dparticle(x0, y0, size, particle+1, orientation_corresp[orientation], axes, representation='arrowed')
        
        
    

    def plot_3Dparticle(self, x0, y0, z0, size, particle, orientation, axes, representation):
        if representation=='arrowed':
            self.plot_3Darrow( x0, y0, z0, size, particle, orientation, axes)
        if representation=="colored":
             self.plot_Cube(x0, y0, z0, size, particle, orientation, axes)


    def plot_particle(self, x0, y0, z0, size, particle, orientation, dimension, axes, representation, order=1):
        if dimension==3 :
            self.plot_3Dparticle(x0, y0, z0, size, particle, orientation, axes, representation)
        if dimension==2:
            if self.lattice.name == 'Hexagonal' :
                """3D hexagons in 2 dimension"""
                self.plot_3Dhexagon(x0, y0, z0, size, particle, orientation, axes, representation=None)
            else :
                self.plot_2Dparticle(x0, y0, size, particle, orientation, axes, representation, order=order)
       

    def plot_one_contact( self, x0, y0, z0, i_edge, radius, structure, axes, LEL, lw=4, emax=None, color_map=None, zorder=None):
            """draw a line of a color corresponding to the energy between particle in x_real, y_real and it's neighbor following edge
            emax is given in absolute value, x0 and y0 are coordinates in real space (not lattice)
            Either LEL or color_map needs to be known to define the colors"""
            
            
            
            if color_map is None:
                e = LEL[structure]
                icol = (e-(-emax))/(2*emax) #number between 0 and 1 that gives the color
                color =  cm.seismic(icol)
            else:
                color = color_map[structure]
                
            
            cos1 = self.cos_vertices[i_edge-1]
            sin1 = self.sin_vertices[i_edge-1]
            cos2 = self.cos_vertices[i_edge]
            sin2 = self.sin_vertices[i_edge]
            
            x1, x2, y1, y2 = x0+cos1*radius, x0+cos2*radius, y0+sin1*radius, y0+sin2*radius
            axes.plot([x1,x2], [y1, y2], color=color, linewidth=lw,  solid_capstyle='round', zorder=zorder )
            # )




"""matrices used for reordering of the lel"""
permut = np.array([[0,1,0,0,0,0],
                   [0,0,1,0,0,0],
                   [0,0,0,1,0,0],
                   [0,0,0,0,1,0],
                   [0,0,0,0,0,1],
                   [1,0,0,0,0,0]])

permut_inv = np.linalg.inv(permut)

anti_diag = np.fliplr(np.eye(6))


class DesignLEL():
    """This class is used to design specific LELs
    The LEL array is built with the energy you choose 
    for a given relative orientation, or pair of particles"""
    
    
    def __init__(self, lattice, n_particles_type):
        
        self.lattice = lattice
        self.n_particles_type = n_particles_type
        self.n_structures = self.lattice.all_N_structures[n_particles_type]
        
        self.N1 = self.lattice.detail_all_N_structures[self.n_particles_type][0]
        
        self.compute_invariance() # see the function for explanation on invariances

    
    
    def cut_LEL(self, LELarg):
        """Takes a LEL that has the dimension of the 1 and 2 particles structures
        that is used with C++ 
        and returns the LEL with only two particles structures"""
        # N1 = self.lattice.detail_all_N_structures[self.n_particles_type][0]
        return(LELarg[self.N1:])
    


   
    def areEquivalent(self, LEL1, LEL2):
        for n_perm in range(self.lattice.n_orientations):
            for d in range(2):
                LELcompare, n, d_chosen = self.reorder_LEL(LEL2, n_permutation=n_perm, direction=d)
                if np.all(LEL1==LELcompare):
                    return(True)
        return(False)

    def LEL_predetermined(self, LELarg):
        """Takes a LEL that has the dimension of the 2 particles structures, 
        and returns the full LEL to write to C++"""
        LEL = np.zeros(self.n_structures)
        # N1 = self.lattice.detail_all_N_structures[self.n_particles_type][0]
        LEL[self.N1:] = LELarg
        return(LEL)
    

    
    
    def interactions_vortex_assembly(self, epsilon_c, epsilon_d, epsilon_f, 
                                     assembly_type = "6-fold"
                                     ):
        
        interaction_map = np.zeros(self.n_structures) 
        
        if assembly_type=="6-fold":
            o1_o2s = [(1,6), (2,3)]
        elif assembly_type=="3-fold":
            o1_o2s = [(2,6), (3,1)]
        elif assembly_type=="fiber":
            o1_o2s = [(1,4),(2,5),(4,1),(3,6)]
        elif assembly_type=="patterned-bulk":
            o1_o2s = [(1,6),(2,3),(4,1),(5,2)]
        elif assembly_type=="4-fold":
            o1_o2s = [(1,2), (1,4), (2,3), (4,3)]
        elif assembly_type=="non-self-limited":
            o1_o2s = [(1,2), (4,3)]
            
        else:
            print("None of the possible assembly_type were chosen")
            print("Choose between 6-fold, 3-fold, fiber, patterned-bulk, 4-fold and non-self-limited")
    
    
        for s in range(self.n_structures):
            state1, state2, edge = self.lattice.structure_to_config_dicts[self.n_particles_type][s]
          
                      
            if (state1!=0 and state2!=0 and edge!=0):
                """Full-full"""
                p1, o1 =  self.lattice.state_to_particle_orientation[state1]
                p2, o2 =  self.lattice.state_to_particle_orientation[state2]
                
                if o1==o2 :
                    """Crystal"""
                    interaction_map[s] = epsilon_c
                
                elif ( (o1==o1_o2s[0][0] and o2==o1_o2s[0][1]) 
                    or (o1==o1_o2s[1][0] and o2==o1_o2s[1][1]) 
                    or (len(o1_o2s)>2 and o1==o1_o2s[2][0] and o2==o1_o2s[2][1]) 
                    or (len(o1_o2s)>3 and o1==o1_o2s[3][0] and o2==o1_o2s[3][1]) ):
                    """Defect line"""
                    interaction_map[s] = epsilon_d
                else :
                    """Forbidden contact"""
                    interaction_map[s] = epsilon_f
        
        return(interaction_map)
    
    
    
    
    
        
    def get_colormap_vortex_assembly(self, color_crystal='mediumaquamarine',
                               color_defect = 'mediumblue', 
                               color_surf = 'tomato', 
                               color_zero = 'white',
                               color_inf = 'darkred', assembly_type="6-fold"):
        
        
        interaction_map = self.interactions_vortex_assembly(-1, -2, 10, assembly_type=assembly_type)[self.N1:]
        s_defect = np.where(interaction_map==-2)[0]
        s_crystal = np.where(interaction_map==-1)[0]
           
        s_sigma =  np.array(range(1,self.N1), dtype=int)
        
        color_map = {}
        color_map[0] = color_zero
        for s in range(1,self.n_structures-self.N1):
            if s in s_defect:
                color_map[s] = color_defect
            elif s in s_crystal :
                color_map[s] = color_crystal
            elif s in s_sigma:
                color_map[s] = color_surf
            else :
                color_map[s] = color_inf
                
        return(color_map)        

        
        
    
    
    
    
    
    
    
    
    def get_x_y_camembert(self, LEL,n_line=6,fiber=False):
        
        if self.n_particles_type==1:
            if not fiber:
                LELex = self.LEL_defect_line(n_line, -1, -2, 0, 10)[self.N1:]
               
            else :
                LELex = self.LEL_fiber_zigzag_2D(-1, -2, 0, 0, 10)[self.N1:]
            s_line = np.where(LELex==-2)[0]
            s_crystal = np.where(LELex==-1)[0]
            # s_line = np.array( [12,15], dtype=int)
            # s_crystal =  np.array([7,14,20], dtype=int)
            
        if self.n_particles_type==2:
            s_line = np.array( [18,21], dtype=int)
            s_crystal = np.array( [13,20,26], dtype=int)
        s_sigma =  np.array(range(1,7), dtype=int)
        
        if  (not np.all(LEL[s_line]==LEL[s_line[0]]) or 
                not np.all(LEL[s_crystal]==LEL[s_crystal[0]]) or
                not np.all(LEL[s_sigma]==LEL[s_sigma[0]])
                
                ):
            print('this LEL does not correspond to camembert')
            return(False)
    
        
        eline = LEL[s_line[0]]
        ecrystal = LEL[s_crystal[0]]
        
        sigma = LEL[s_sigma[0]]
        x = eline-2*sigma
        y = ecrystal - 2*sigma
        return(x,y, 6*sigma)



    
    def compute_invariance(self):
        """this function creates the invariances table
        invariances[i] is a vector such that the scalar product c.invariances[i] is always a constant
        c.invariances[0] should be one (conservation of the bonds)
        c.inverainces[i!=0] shoud be density/3/nparticles_type
           this is the conservation of the number of faces"""        
        
        
        invariances = np.zeros((self.lattice.n_neighbors*self.n_particles_type+1, self.n_structures-self.N1), dtype=int)
        
        for s in range(self.N1, self.n_structures):
             p1, f1, p2, f2 = self.lattice.get_faces_in_contact(s, self.n_particles_type)             
             """the faces in contact :  face 3 of particle 2 is labeled 2*n_faces + 3"""
             f1_tot = p1*self.lattice.n_neighbors + f1
             f2_tot = p2*self.lattice.n_neighbors + f2
             if not (f1_tot==0 and f1_tot==0):
                     invariances[f1_tot,s-self.N1]+=1
                     invariances[f2_tot,s-self.N1]+=1
        invariances[0,:]=np.ones(self.n_structures-self.N1)
        self.invariances = invariances
    
    
    
    
    def transform_in_matrix(self, LELshort, p1=0, p2=0, symmetric=True, emptyFullAlso=False):
        """Takes the short verision of LEL, without N1P structures"""
        LEL = self.LEL_predetermined(LELshort)
        structures_matrix = self.lattice.get_matrix_structures(p1=p1, p2=p2, 
                                                               symmetric=symmetric, 
                                                               n_type=self.n_particles_type,
                                                               emptyFullAlso=emptyFullAlso)
        n, n = structures_matrix.shape
        matrix = np.zeros((n,n))
        for i in range(n):
            for j in range (n):
                matrix[i,j] = LEL[structures_matrix[i,j]]
                
        return(matrix)
    
    
    def transform_from_matrix(self, smallMatrix, p1=0, p2=0, symmetric=True):
        """Takes the nface*nface matrix and returns LEL"""
        
        LEL = np.zeros(self.n_structures)
        structures_matrix = self.lattice.get_matrix_structures(p1=p1, p2=p2, symmetric=symmetric,n_type=self.n_particles_type)
        n, n = structures_matrix.shape
        for i in range(n):
            for j in range (n):
                LEL[structures_matrix[i,j]] = smallMatrix[i,j] 
        return(self.cut_LEL(LEL))
        
    
    
    
    def reorder_LEL(self, LELshort, n_permutation=None, direction=None, symmetric=True):
        """takes short LEL and return short LEL
        if n_permutation and direction are given, the reordering is predetermined"""
        matrix = self.transform_in_matrix(LELshort, symmetric=symmetric)
        empty_fulls = LELshort[1:self.lattice.n_neighbors+1]
        
        """Permutation of the matrix"""
        matrix = matrix.copy()
        empty_fulls.copy()
        
        if n_permutation is None :
            dim = matrix.shape[0]
            eig, eigv = np.linalg.eigh(matrix)
            eigv = eigv.T
            order = np.argsort(eig)
            # print('order',order)
            # print('eig', eig[order])
            eigv = eigv[order]
            eig = eig[order]
            
            target_best_i = 0
            
            """choosing the most negative lambda*u_i"""
            # mult = np.zeros((dim,dim))
            # for i in range(dim):
            #     mult[i] = eigv[i]*eig[i] #this is a vector
            # min_mult = np.argmin(mult)
            # absolute_max, best_i = np.unravel_index(min_mult, (dim,dim))
                
            # all_mult = np.unique(mult)
            # # print('min_mult and next', all_mult[0],all_mult[1])
            # abs_max_eigenvector = np.abs(eigv[absolute_max])
            # direction = np.argmin([abs_max_eigenvector[(best_i-1)%dim], abs_max_eigenvector[(best_i+1)%dim]]) #this is 0 or 1
           
            
            """choosing biggest projection on (1,0,0..) for the biggest eigenvector"""
            absolute_max = np.argmin(np.abs(eig))
            # eigenvector of the lowest eigenvalue :
            abs_max_eigenvector = np.abs(eigv[absolute_max])
            # print('eigv', eigv[absolute_max])
            # print('eigvAbs', abs_max_eigenvector)
            
            #We choose the permutation of the eigenvector 
            # with the biggest projection on [1,0,0,...] 
            best_i = np.argmin(abs_max_eigenvector) 
            # print('absolute_max', 'best_i',absolute_max,  best_i)
            
            # We choose the permutation of the eigenvector 
            # with the biggest projection on [0,1,0,...] 
            # direction is 0 or 1    
            direction = np.argmin([abs_max_eigenvector[(best_i-1)%dim], abs_max_eigenvector[(best_i+1)%dim]]) #this is 0 or 1
           
            
            """Choosing with the faces"""
            # face_stickiness = np.average(matrix, axis=1)
            # best_i = np.argmin(face_stickiness)
            
            # #direction is 0 or 1           
            # direction = np.argmin([face_stickiness[(best_i-1)%dim], face_stickiness[(best_i+1)%dim]]) 
            
            
            # print('best_i', best_i)
            # print('direction', direction)
            if direction==0 :
                n_permutation = np.abs(target_best_i - best_i)
            else :
                n_permutation = (dim-1) - np.abs(target_best_i - best_i)
                
    
        # we will do n_permutation cyclic permutation
        M = np.linalg.matrix_power(permut, n_permutation)
        Minv = np.linalg.inv(M)
        new_matrix = M.dot(matrix.dot(Minv))
        
        new_empty_fulls = M.dot(empty_fulls)
        #if needed, we also do the miror permutation
        if direction == 1 :
            new_matrix = anti_diag.dot(new_matrix.dot(anti_diag))
            new_empty_fulls = anti_diag.dot(empty_fulls)
        
        
        newLEL = np.zeros(len(LELshort))
        
        newLELfullfull =  self.transform_from_matrix(new_matrix, symmetric=symmetric)
        newLEL =newLELfullfull
        newLEL[0] = LELshort[0]
        newLEL[1:self.lattice.n_neighbors+1] = new_empty_fulls
        
        return(newLEL, n_permutation, direction)
    
    

class LELrepresentation():
    
    def __init__(self, lattice, n_particles_type):
        self.lattice = lattice
        self.n_particles_type = n_particles_type
        self.n_structures = self.lattice.all_N_structures[n_particles_type]
        
        self.designs = DesignLEL(self.lattice, self.n_particles_type)
        
        
    def get_colormatrixLEL(self, LEL, p1=0, p2=0,colormap=cm.seismic, 
                           colordict=None, maxi=None,
                           structures_in_grey=None,
                           grey_color='gray',
                           symmetric=True,
                           emptyFullAlso = False):
        """Takes the short verision of LEL, without N1P structures"""

        
        matrix = self.designs.transform_in_matrix(LEL, p1=p1, p2=p2, symmetric=symmetric, emptyFullAlso=emptyFullAlso)
        n = self.lattice.n_orientations + emptyFullAlso
        N1 = self.lattice.detail_all_N_structures[self.n_particles_type][0]
        
        
        if structures_in_grey is not None:
            if maxi is None and LEL is not None:
                maxi = max(-np.min(LEL), np.max(LEL))   
            structures_matrix = self.lattice.get_matrix_structures(p1=p1, p2=p2, 
                                                                   symmetric=symmetric,
                                                                   n_type=self.n_particles_type,
                                                                   emptyFullAlso=emptyFullAlso)-N1
            colorMatrix = np.zeros((n,n,4))
            for i in range (n):
                for j in range (n):
                    s = structures_matrix[i,j]
                    if s in structures_in_grey:
                       colorMatrix[i,j]  =  matplotlib.colors.to_rgba(grey_color)
                    else :
                        if LEL is None :
                            colorMatrix[i,j]  =  matplotlib.colors.to_rgba('white')
                        else :
                            J = matrix[i][j]
                            level = 0.5 + J/maxi
                            colorMatrix[i,j] = colormap(level)
                            
        
        elif colordict is not None :
            structures_matrix = self.lattice.get_matrix_structures(p1=p1, p2=p2, 
                                                                   symmetric=symmetric,
                                                                   n_type=self.n_particles_type,
                                                                   emptyFullAlso=emptyFullAlso)-N1
            colorMatrix = np.zeros((n,n,4))
            for i in range (n):
                for j in range (n):
                    colorMatrix[i,j] = matplotlib.colors.to_rgba(colordict[structures_matrix[i,j]])
        else :
            if maxi is None :
                maxi = max(-np.min(LEL), np.max(LEL))
            colorMatrix = np.zeros((n,n,4))
            for i in range (n):
                for j in range (n):
                    J = matrix[i][j]
                    level = 0.5 + J/maxi
                    colorMatrix[i,j] = colormap(level)
                
        return(colorMatrix)
    
    def get_colormatrixDensity(self, vector, p1=0, p2=0,
                               colormap=cm.seismic,maxi=None,
                               symmetric=True):
        """Takes the short verision of c, without N1P structures"""

        
        matrix = self.designs.transform_in_matrix(vector, p1=p1, p2=p2, symmetric=symmetric).copy()
        n = self.lattice.n_orientations
        N1 = self.lattice.detail_all_N_structures[self.n_particles_type][0]
        
        if maxi is None :
            maxi = np.max(vector[N1:])
        # print(maxi)
        colorMatrix = np.zeros((n,n,4))
        for i in range (n):
            for j in range (n):
                J = matrix[i][j]
                if i!=j:
                    J = J/2
                level = 0.5 - J/maxi
                colorMatrix[i,j] = colormap(level)

        return(colorMatrix)
    
    def plot_matrix(self, vector, vectorIsLEL=True, 
                    p1=0, p2=0,
                    axes=None, 
                    x0=0, y0=0, dx=1,dy=1,
                    fullMatrix = True, #if False, plot only the triangular inferior part
                    emptyFullAlso = False, #if true, plot also the structures empty-full
                    lw = 0.5,axes_off=True,
                    colormap=cm.seismic, 
                    colordict=None, #this is to colorcode the structures, not the energy level
                    structures_in_grey = None,
                    grey_color = 'gray',
                    arangeAxis = True,
                    maxi=None, 
                    margin=2, symmetric=True, 
                    plot_index_particle=False):
        """Takes the short verision of LEL, without N1P structures"""
        
        if vectorIsLEL :
            colorMatrix = self.get_colormatrixLEL(vector, p1=p1, p2=p2,
                                                  colormap=colormap, colordict=colordict, 
                                                  structures_in_grey=structures_in_grey,
                                                  grey_color=grey_color,
                                                  maxi=maxi, symmetric=symmetric, emptyFullAlso=emptyFullAlso)
        else  :
            colorMatrix = self.get_colormatrixDensity(vector, p1=p1, p2=p2,
                                                      colormap=colormap,maxi=maxi
                                                      , symmetric=symmetric)
                
        n,n, colordimension=colorMatrix.shape
    
        if axes is None :
            f,axes = plt.subplots()
    
        for i in range(n):
            if fullMatrix :
                jlist = range(0,n) 
            else :
                jlist = range(i,n)
            for j in jlist:
                color = colorMatrix[i,j]
                square = mpatches.RegularPolygon((x0+i*dx,y0-j*dy), 4, orientation=np.pi/4, 
                                                 radius=dx/np.sqrt(2),  color=color,ec='black', 
                                                 lw=lw)
                axes.add_artist(square)
    
        if arangeAxis :
            margin = margin
            axes.set_xlim(x0-margin, x0+(n-1)*dx+margin)
            axes.set_ylim(y0-(n-1)*dx - margin, y0+margin)
            axes.set_aspect('equal', 'box')  
            axes.set_xticks([])
            axes.set_yticks([])#,[])
            
        
        else:
            n1, n2 = -2, 7
            axes.set_xlim(x0+n1*dx, x0+n2*dx)
            axes.set_ylim(y0-n2*dx, y0-n1*dx)
            axes.set_aspect('equal', 'box')
        if axes_off :
            axes.axis('off')
            
        """plot matrix indices as particles"""
        if plot_index_particle:
            
            PR = ParticleRepresentation(self.lattice)
            
            dxP = dx
            rx = 0.7*dx
            ry=rx
            
            if emptyFullAlso:
                PR.plot_2Darrowed_particle(x0,y0+dxP, rx, 0, 0, axes, edgecolor='black')
                PR.plot_2Darrowed_particle(x0-dxP, y0, rx, 0, 0, axes, edgecolor='black')
            # else:
            #     x0 = x0-dxP
            #     y0 = y0-dx
            
            # plt.scatter(x0,y0+dxP, color='gray')
            if self.lattice.name=='Triangular':
                o1s = [0,1,2,3,4,5,6]
                o2s = [0,4,5,6,1,2,3]
                my_range = range(1,7)
            if self.lattice.name=='Square':
                o1s = [0,1,2,3,4]
                o2s = [0,3,4,1,2,]
                my_range = range(1,5)
                
            for i in my_range:
                o = o2s[i]
                PR.plot_2Darrowed_particle(x0=x0+(i-1*(not emptyFullAlso))*dx, y0=y0+dxP, size=rx, 
                                           orientation=o,particle=0, axes=axes, edgecolor='black')
            for j in my_range:
                o = o1s[j]
                PR.plot_2Darrowed_particle(x0=x0-dxP, y0=y0-(j-1*(not emptyFullAlso))*dx, size=rx, 
                                           orientation=o, particle=0, axes=axes, edgecolor='black')
    
    
            

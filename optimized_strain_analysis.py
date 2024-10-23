from ase import Atoms
import ase
from ase.geometry.analysis import Analysis
import matplotlib.pyplot as plt
import numpy as np
from  matplotlib.colors import LinearSegmentedColormap
from  matplotlib.colors import Normalize
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.cm import ScalarMappable
from ase import build
import matplotlib.ticker as ticker
from matplotlib.collections import LineCollection

class StrainAnalysis:
    def __init__(self):
        self.bond_z_change = None
        self.bondlist_left_coor_x = None
        self.bondlist_right_coor_x = None
        self.bondlist_left_coor_y = None
        self.bondlist_right_coor_y = None
        self.str_opt_x = None
        self.str_opt_y = None
        self.bond_diff = None
        self.z_change = None
    def d_interlayer(self,layer1,layer2):
        coor_1 = [i[2] for i in layer1.get_positions()]
        coor_2 = [i[2] for i in layer2.get_positions()]
        interlayer_d = abs(np.mean(coor_1)-np.mean(coor_2))
        return interlayer_d
    def d_vacuum(self,layer):
        coor = [i[2] for i in layer.get_positions()]
        vacuum = abs(np.mean(coor))
        return vacuum
    
    def layer_sep(self,geometry,peoridicity): #geometry is Atoms object from ase
        bond_list = Analysis(geometry).all_bonds[0]
        collection = []
        single_layer = []
        trig = bond_list[0]
        duplicate_check = []
        single_layer = trig
        while True:
            trig_copy = trig
            trig = []
            for j in trig_copy:
                for m in bond_list[j]:
                    trig.append(m)
            trig = list(set(trig))
            for atom in trig:
                if atom in single_layer:
                    duplicate_check.append(0)
                else:
                    single_layer.append(atom)
            if len(duplicate_check) == len(trig):
                output = list(set(single_layer))
                duplicate_check = []
                for ind in output:
                    bond_list[ind] = 0
                single_layer = []
                for i in bond_list:
                    if i!=0:
                        trig = i
                for i,item in enumerate(output):
                    output[i] = geometry[item]
                if peoridicity == True:
                    output_atoms = Atoms(output,pbc = True)
                    output_atoms.set_cell(geometry.get_cell())
                    collection.append(output_atoms)
                else:
                    output_atoms = Atoms(output)
                    collection.append(output_atoms)
            else:
                duplicate_check = []
            check = np.array(bond_list,dtype=object)
            if (check== 0).all() == True:
                break
        return collection
    
    def flake_generator(self,geometry,x,y,z):
        #Atoms.center(geometry)
        P = np.zeros([3,3])
        P[0][0] = x
        P[1][1] = y
        P[2][2] = z
        supercell = ase.build.make_supercell(geometry,P)
        supercell.set_pbc(False)
        supercell.set_cell(None)
        return supercell
    
    def mirror(self,geometry,plane):
        positions = geometry.get_positions()
        x = list(map(lambda x:x[0],positions))
        y = list(map(lambda y:y[1],positions))
        z = list(map(lambda z:z[2],positions))
        x_mean = np.mean(x)
        y_mean = np.mean(y)
        z_mean = np.mean(z)
        if plane == 'xy':
            new_positions = list(map(lambda coor:[coor[0],coor[1],2*z_mean-coor[2]],positions))
            geometry.positions = new_positions
        if plane == 'xz':
            new_positions = list(map(lambda coor:[coor[0],2*coor[1]-y_mean,coor[2]],positions))
            geometry.positions = new_positions
        if plane == 'yz':
            new_positions = list(map(lambda coor:[2*x_mean-coor[0],coor[1],coor[2]],positions))
            geometry.positions = new_positions
        return 
    
    def prepare_G(self,str_orig,str_opt,refernece_d):
#python simple plot for strain
#prepare data of structure
#structures must be flake
#str_orig,str_opt are atoms objects of flakes
        str_orig_x = []
        str_orig_y = []
        str_orig_z = []
        str_opt_x = []
        str_opt_y = []
        str_opt_z = []
        z_change = []
        bond_length_orig = []
        bond_length_opt = []
        z_atoms_up = []
        z_atoms_down = []
        z_atoms_neutral = []
        bondlist_left = []
        bondlist_right = []
        bondlist_left_coor_x = []
        bondlist_right_coor_x = []
        bondlist_left_coor_y = []
        bondlist_right_coor_y = []
        bondlist_left_coor_z = []
        bondlist_right_coor_z = []
        bond_z_change = []

        str_orig_coor = str_orig.get_positions()
        str_opt_coor = str_opt.get_positions()
        for c1,c2 in zip(str_orig_coor,str_opt_coor):
            str_orig_x.append(c1[0])
            str_orig_y.append(c1[1])
            str_orig_z.append(c1[2])
            str_opt_x.append(c2[0])
            str_opt_y.append(c2[1])
            str_opt_z.append(c2[2])
            
        str_orig_bond = Analysis(str_orig).all_bonds
        str_opt_bond = Analysis(str_opt).all_bonds
        
#get the change along z axis and within horizontal plane        
        for i in str_opt_z:
            z_change.append(i-refernece_d)
            #z_change.append(i-refernece_d-np.average(str_opt_z))
           
        for count,coor in enumerate(z_change):
            if coor > 0:
                z_atoms_up.append(count)
            elif coor < 0:
                z_atoms_down.append(count)
            elif coor == 0:
                z_atoms_neutral.append(count)
    
#bond length change
        for atom,bonded_array in zip(range(len(str_orig_z)),str_orig_bond[0]):
            for bonded in bonded_array:
                bond_length_orig.append(str_orig.get_distance(atom,bonded))
        for atom,bonded_array in zip(range(len(str_opt_z)),str_opt_bond[0]):
            for bonded in bonded_array:
                bond_length_opt.append(str_opt.get_distance(atom,bonded))

        bond_diff = [i-j for i,j in zip(bond_length_opt,bond_length_orig)]


##create list for bonds of graphene, compare the center with average to determine up or down        

        for i in range(len(str_opt_bond[0])):
            for index in str_opt_bond[0][i]:
                bondlist_left.append(i)
                bondlist_right.append(index)
        for i,j in zip(bondlist_left,bondlist_right):
            bondlist_left_coor_x.append(str_opt_x[i])
            bondlist_right_coor_x.append(str_opt_x[j])
            bondlist_left_coor_y.append(str_opt_y[i])
            bondlist_right_coor_y.append(str_opt_y[j])
            bondlist_left_coor_z.append(str_opt_z[i])
            bondlist_right_coor_z.append(str_opt_z[j])
        for i,j in zip(bondlist_left_coor_z,bondlist_right_coor_z):
            bond_z_change.append((i+j)/2-refernece_d)
            #bond_z_change.append((i+j)/2-np.average(str_opt_z))
            
        self.bond_z_change = bond_z_change
        self.bondlist_left_coor_x = bondlist_left_coor_x
        self.bondlist_right_coor_x = bondlist_right_coor_x
        self.bondlist_left_coor_y = bondlist_left_coor_y
        self.bondlist_right_coor_y = bondlist_right_coor_y
        self.str_opt_x = str_opt_x
        self.str_opt_y = str_opt_y
        self.bond_diff = bond_diff
        self.z_change = z_change
        return 
          
    def prepare_COF(self,str_orig,str_opt):
#python simple plot for strain
#prepare data of structure
#structures must be flake
#str_orig,str_opt are atoms objects of flakes
        str_orig_x = []
        str_orig_y = []
        str_orig_z = []
        str_opt_x = []
        str_opt_y = []
        str_opt_z = []
        z_change = []
        bond_length_orig = []
        bond_length_opt = []
        z_atoms_up = []
        z_atoms_down = []
        z_atoms_neutral = []
        bondlist_left = []
        bondlist_right = []
        bondlist_left_coor_x = []
        bondlist_right_coor_x = []
        bondlist_left_coor_y = []
        bondlist_right_coor_y = []
        bondlist_left_coor_z = []
        bondlist_right_coor_z = []
        bond_z_change = []

        str_orig_coor = str_orig.get_positions()
        str_opt_coor = str_opt.get_positions()
        for c1,c2 in zip(str_orig_coor,str_opt_coor):
            str_orig_x.append(c1[0])
            str_orig_y.append(c1[1])
            str_orig_z.append(c1[2])
            str_opt_x.append(c2[0])
            str_opt_y.append(c2[1])
            str_opt_z.append(c2[2])
            
        str_orig_bond = Analysis(str_orig).all_bonds
        str_opt_bond = Analysis(str_opt).all_bonds
        
#get the change along z axis and within horizontal plane        
        for i in str_opt_z:
            z_change.append(i-np.average(str_opt_z))

        for count,coor in enumerate(z_change):
            if coor > 0:
                z_atoms_up.append(count)
            elif coor < 0:
                z_atoms_down.append(count)
            elif coor == 0:
                z_atoms_neutral.append(count)
    
#bond length change
        for atom,bonded_array in zip(range(len(str_orig_z)),str_orig_bond[0]):
            for bonded in bonded_array:
                bond_length_orig.append(str_orig.get_distance(atom,bonded))
        for atom,bonded_array in zip(range(len(str_opt_z)),str_opt_bond[0]):
            for bonded in bonded_array:
                bond_length_opt.append(str_opt.get_distance(atom,bonded))

        bond_diff = [i-j for i,j in zip(bond_length_opt,bond_length_orig)]


##create list for bonds of graphene, compare the center with average to determine up or down        

        for i in range(len(str_opt_bond[0])):
            for index in str_opt_bond[0][i]:
                bondlist_left.append(i)
                bondlist_right.append(index)
        for i,j in zip(bondlist_left,bondlist_right):
            bondlist_left_coor_x.append(str_opt_x[i])
            bondlist_right_coor_x.append(str_opt_x[j])
            bondlist_left_coor_y.append(str_opt_y[i])
            bondlist_right_coor_y.append(str_opt_y[j])
            bondlist_left_coor_z.append(str_opt_z[i])
            bondlist_right_coor_z.append(str_opt_z[j])
        for i,j in zip(bondlist_left_coor_z,bondlist_right_coor_z):
            bond_z_change.append((i+j)/2-np.average(str_opt_z))
            
        self.bond_z_change_cof = bond_z_change
        self.bondlist_left_coor_x_cof = bondlist_left_coor_x
        self.bondlist_right_coor_x_cof = bondlist_right_coor_x
        self.bondlist_left_coor_y_cof = bondlist_left_coor_y
        self.bondlist_right_coor_y_cof = bondlist_right_coor_y
        self.str_opt_x_cof = str_opt_x
        self.str_opt_y_cof = str_opt_y
        self.bond_diff_cof = bond_diff
        self.z_change_cof = z_change    
        return
    
    def out_plane(self, color, transparency=0.7, linewidth=0.5, vmin=None, vmax=None, labelsize=14, fontsize=16, colorbar_tick_locator=0.1):
        # Convert bond coordinates to segments format for LineCollection
        segments = np.array([
            [[x1, y1], [x2, y2]] for x1, y1, x2, y2 in zip(
                self.bondlist_left_coor_x,
                self.bondlist_left_coor_y,
                self.bondlist_right_coor_x,
                self.bondlist_right_coor_y
            )
        ])
        
        # Set up color mapping
        if vmin is None:
            vmin = min(self.bond_z_change)
        if vmax is None:
            vmax = max(self.bond_z_change)
            
        cmap = LinearSegmentedColormap.from_list('rgba', color, N=256)
        norm = Normalize(vmin=vmin, vmax=vmax)
        
        # Create figure and axis
        fig, ax = plt.subplots()
        ax.set_aspect('equal', adjustable='box')
        
        # Create LineCollection for bonds
        lc = LineCollection(segments, linewidths=linewidth, colors=cmap(norm(self.bond_z_change)))
        ax.add_collection(lc)
        
        # Plot atoms using scatter (more efficient than multiple plot calls)
        scatter = ax.scatter(self.str_opt_x, self.str_opt_y, s=1, 
                           c=self.z_change, cmap=cmap, alpha=transparency,
                           vmin=vmin, vmax=vmax, zorder=5)
        
        # Set plot limits based on data
        ax.set_xlim(np.min(self.str_opt_x) - 1, np.max(self.str_opt_x) + 1)
        ax.set_ylim(np.min(self.str_opt_y) - 1, np.max(self.str_opt_y) + 1)
        
        # Add colorbar
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad="3%")
        cb = plt.colorbar(scatter, cax=cax)
        cb.locator = ticker.MultipleLocator(base=colorbar_tick_locator)
        cb.update_ticks()
        cb.ax.tick_params(labelsize=labelsize)
        cb.set_label('\u0394'+'d'+' / '+'$\AA$', fontsize=fontsize)
        cb.ax.yaxis.set_ticks_position('right')
        
        ax.set_xticks([])
        ax.set_yticks([])
        
        return fig, ax
    
    def in_plane(self, color, linewidth, vmin=None, vmax=None):
        # Convert bond coordinates to segments format
        segments = np.array([
            [[x1, y1], [x2, y2]] for x1, y1, x2, y2 in zip(
                self.bondlist_left_coor_x,
                self.bondlist_left_coor_y,
                self.bondlist_right_coor_x,
                self.bondlist_right_coor_y
            )
        ])
        
        # Set up color mapping
        if vmin is None:
            vmin = min(self.bond_diff)
        if vmax is None:
            vmax = max(self.bond_diff)
            
        cmap_bond = LinearSegmentedColormap.from_list('rgba', color, N=256)
        norm = Normalize(vmin=vmin, vmax=vmax)
        
        # Create figure and axis
        fig, ax = plt.subplots()
        ax.set_aspect('equal', adjustable='box')
        
        # Create and add LineCollection
        lc = LineCollection(segments, linewidths=linewidth, cmap=cmap_bond, norm=norm)
        lc.set_array(np.array(self.bond_diff))
        ax.add_collection(lc)
        
        # Set plot limits
        ax.set_xlim(np.min(self.bondlist_left_coor_x) - 1, np.max(self.bondlist_right_coor_x) + 1)
        ax.set_ylim(np.min(self.bondlist_left_coor_y) - 1, np.max(self.bondlist_right_coor_y) + 1)
        
        # Add colorbar
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad="3%")
        cb = plt.colorbar(lc, cax=cax)
        cb.ax.tick_params(labelsize=12)
        cb.set_label('\u0394'+'d'+' / '+'$\AA$', fontsize=13)
        
        ax.set_xticks([])
        ax.set_yticks([])
        
        return fig, ax

#plot COF frame on top of graphene
    def add_COF(self,fig,ax,color,name,path,show_atoms = 'False',clear_after = 'True',transparancy = 0.7,zorder = 10,linewidth = 0.5):
        ax2 = fig.add_subplot(111,sharex = ax, sharey = ax)
        ax2.set_position(ax.get_position())
        ax2.patch.set_visible(False)
        for i in range(len(self.bondlist_left_coor_x_cof)):
            ax2.plot([self.bondlist_left_coor_x_cof[i],self.bondlist_right_coor_x_cof[i]],[self.bondlist_left_coor_y_cof[i],self.bondlist_right_coor_y_cof[i]],linewidth=linewidth,c=color,alpha = transparancy,zorder = zorder)
        if show_atoms == 'True':
            ax2.scatter(self.str_opt_x_cof, self.str_opt_y_cof,s=1,c=color,alpha = transparancy,zorder = zorder)
        fig.savefig(f'{path}{name}.png',dpi = 600, bbox_inches = 'tight')
        if clear_after == 'True':
            ax2.cla()
        return
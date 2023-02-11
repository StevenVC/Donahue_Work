"""
Created on Mon Jul 11 12:43:11 2022

@author: steven
"""

# import sherpa
import sherpa.astro.ui as sh
import ciao_contrib.all

def gAbsorber_powerlaw(id, index, nH):
    '''
    Identifier: 'gp'
    
    Souce model with components:
        galactic absorber (*)
        powerlaw (+)
        
    nH: galactic absorbption (float) [10^22 atoms / cm^2]
    '''
    i_n = str(index)

    # set source model
    model = f'xsphabs.abs{i_n} * powlaw1d.p{i_n}'
    sh.set_source(id, model)
    
    # set model components
    sh.set_par(f'abs{i_n}.nH', val=nH, frozen=True)

    return


def gAbsorber_zAbsorber_powerlaw(id, index, nH, redshift):
    '''
    Description:
        This function sets the source model in the sherpa environment with the following properties
        
        Identifier: 'gzp'
        
        Source model with components:
            galactic absorber (*)
            redshifted absorber (*)
            powerlaw (+)
        
    Inputs:
        nH (float): 
            galactic absorbption; units [10^22 atoms / cm^2]
        
        redshift (float): 
            redshift of zAbsorber; units [z]
            
    Returns:
        None
    '''
    
    # set source model
    i_n = str(index)

    model = f'xsphabs.abs{i_n} * xszphabs.zabs{i_n} * powlaw1d.p{i_n}'
    sh.set_source(id, model)
    
    # set model components
    sh.set_par(f'abs{i_n}.nH', val=nH, frozen=True)

    sh.set_par(f'zabs{i_n}.nH', val=0, min=0, frozen=False)
    sh.set_par(f'zabs{i_n}.redshift', val=redshift, frozen=False)


def gAbsorber_zAbsorber_emissionSpec_powerlaw(id, index, nH, redshift, kT):
    '''
    Description:
        This function sets the source model in the sherpa environment with the following properties
        
        Identifier: 'gzep'
        
        Source model with components:
            galactic absorber (*)
            redshifted absorber (*)
            APEC emission spectrum (+)
            powerlaw (+)
        
    Inputs:
        nH (float): 
            galactic absorbption; units [10^22 atoms / cm^2]
        
        redshift (float): 
            redshift of zAbsorber; units [z]
            
        kT (float): 
            temperature of emission spectrum; units [keV]
            
    Returns:
        None
        
    '''
    
    # set source model
    i_n = str(index)

    model = f'xsphabs.abs{i_n} * xszphabs.zabs{i_n} * (xsapec.apec{i_n} + powlaw1d.p{i_n})'
    sh.set_source(id, model)

    # set model components
    sh.set_par(f'abs{i_n}.nH', val=nH, frozen=True)

    sh.set_par(f'zabs{i_n}.nH', val=0, min=0, frozen=False)
    sh.set_par(f'zabs{i_n}.redshift', val=redshift, frozen=False)

    sh.set_par(f'apec{i_n}.kT', val=kT, frozen=True)
    sh.set_par(f'apec{i_n}.Abundanc', val=1, frozen=True)
    sh.set_par(f'apec{i_n}.redshift', val=redshift, frozen=True)
    sh.set_par(f'apec{i_n}.norm', frozen=False)


def gAbsorber_pAbsorber_powerlaw(id, index, nH, redshift, CvrFract):
    '''
    Description:
        This function sets the source model in the sherpa environment with the following properties
        
        Identifier: 'gpp'
        
        Source model with components:
            galactic absorber (*)
            partial covering absorber (*)
            powerlaw (+)
        
    Inputs:
        nH (float): 
            galactic absorbption; units [10^22 atoms / cm^2]
        
        redshift (float): 
            redshift of zAbsorber; units [z]
            
        CvrFract (float):
            partial covering fraction; units []
    
    Returns:
        None
    '''
    
    # set source model
    i_n = str(index)
    
    model = f'xsphabs.abs{i_n} * xszpcfabs.pabs{i_n} * powlaw1d.p{i_n}'
    
    sh.set_source(id, model)
    
    # set model components
    sh.set_par(f'abs{i_n}.nH', val=nH, frozen=True)

    sh.set_par(f'pabs{i_n}.CvrFract', val=CvrFract, frozen=False)
    sh.set_par(f'pabs{i_n}.redshift', val=redshift, frozen=True)
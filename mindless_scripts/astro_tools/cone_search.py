import click
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord, Distance
import astropy.units as u
from astropy.time import Time
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
from datetime import datetime
import sys
import os
import warnings

# Add the project root to sys.path to allow imports of mindless_scripts package
# This is necessary when running the script directly from the file system
try:
    # Try importing normally first (in case installed or in path)
    from mindless_scripts.astro_tools.unicoord import unicoord
    from mindless_scripts.astro_tools.prop_motion import prop_motion
except ImportError:
    # If that fails, try adding the project root to sys.path
    current_dir = os.path.dirname(os.path.abspath(__file__))
    # Go up two levels: astro_tools -> mindless_scripts (package) -> mindless_scripts (root)
    # Wait, structure is ROOT/mindless_scripts/astro_tools/cone_search.py
    # We want ROOT in path.
    # dirname -> astro_tools
    # dirname(dirname) -> mindless_scripts (package)
    # dirname(dirname(dirname)) -> ROOT
    
    project_root = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    if project_root not in sys.path:
        sys.path.insert(0, project_root)

    try:
        from mindless_scripts.astro_tools.unicoord import unicoord
        from mindless_scripts.astro_tools.prop_motion import prop_motion
    except ImportError:
        # Fallback for local imports if structure is different or running from same dir
        # prop_motion.py uses absolute import 'from mindless_scripts...', so we MUST have ROOT in path for it to work.
        # If we are here, it means even with ROOT in path it failed, or we couldn't find ROOT.
        # But let's try local import as last resort for unicoord, but prop_motion might still fail if it doesn't use relative imports.
        # prop_motion.py content: "from mindless_scripts.astro_tools.unicoord import unicoord"
        # So we absolutely need 'mindless_scripts' to be resolvable.
        sys.path.append(current_dir)
        from unicoord import unicoord
        from prop_motion import prop_motion

# Radio catalogues from the reference notebook
RADIO_CATALOGUES = {
    'nvss': ['VIII/65', 'NVSS'],
    'sumss': ['VIII/81B', 'SUMSS'],
    'wendker_2001': ['VIII/99', 'Wendker'],
    'srsc': ['J/other/PASA/41.84', 'SRSC'],
    'barrett_mcv': ['J/AJ/154/252', 'Barrett'],
    'gleamx': ['VIII/113', 'GLEAM-X'],
    'vlass': ['J/ApJS/255/30', 'VLASS'],
    'morx': ['V/158', 'MORX'],
    'first': ['VIII/92', 'FIRST'],
    'atca_smc': ['J/MNRAS/335/1085', 'ATCA_SMC']
}

def perform_radio_search(coord, radius_arcsec):
    """
    Perform a cone search across defined radio catalogues.
    """
    v = Vizier(columns=['*', '+_r'])
    results = []
    
    for cat_key, (cat_id, cat_name) in RADIO_CATALOGUES.items():
        try:
            res = v.query_region(coord, radius=radius_arcsec * u.arcsec, catalog=cat_id)
            if res and len(res) > 0:
                df = res[0].to_pandas()
                
                for _, row in df.iterrows():
                    ra_col = [c for c in row.index if 'RA' in c.upper() and 'J2000' in c.upper()]
                    dec_col = [c for c in row.index if 'DE' in c.upper() and 'J2000' in c.upper()]
                    
                    ra = row[ra_col[0]] if ra_col else np.nan
                    dec = row[dec_col[0]] if dec_col else np.nan
                    
                    if pd.isna(ra) or pd.isna(dec):
                        continue

                    dist = row['_r'] if '_r' in row else np.nan
                    
                    source_id = f"{cat_name} source"
                    
                    results.append({
                        'ID': source_id,
                        'RA': ra,
                        'Dec': dec,
                        'Distance_arcsec': dist * 60 if dist < 1 else dist, # Assume deg if < 1 (unlikely for _r in Vizier usually arcsec/deg mixed but let's trust _r)
                        # Actually Vizier _r is usually in degrees if not specified? 
                        # Let's recalculate distance to be safe.
                        'Catalogue': cat_name,
                        'PMRA': np.nan,
                        'PMDec': np.nan,
                        'Plx': np.nan
                    })
        except Exception as e:
            pass

    # Recalculate distances for consistency
    final_results = []
    for res in results:
        match_coord = SkyCoord(ra=res['RA'], dec=res['Dec'], unit=(u.deg, u.deg), frame='icrs')
        res['Distance_arcsec'] = coord.separation(match_coord).arcsec
        final_results.append(res)

    return pd.DataFrame(final_results)

def perform_simbad_search(coord, radius_arcsec):
    """
    Perform a cone search using Simbad.
    """
    custom_simbad = Simbad()
    custom_simbad.add_votable_fields('ra(d)', 'dec(d)', 'pmra', 'pmdec', 'plx', 'id(1)')
    
    try:
        res = custom_simbad.query_region(coord, radius=radius_arcsec * u.arcsec)
        if res is None:
            return pd.DataFrame()
        
        df = res.to_pandas()
        results = []
        for _, row in df.iterrows():
            # Handle column name changes in astroquery/simbad
            # Try new names first (ra, dec), then old names (RA_d, DEC_d), then uppercase (RA, DEC)
            ra_val = row.get('ra', row.get('RA', row.get('RA_d', np.nan)))
            dec_val = row.get('dec', row.get('DEC', row.get('DEC_d', np.nan)))
            
            match_coord = SkyCoord(ra=ra_val, dec=dec_val, unit=(u.deg, u.deg), frame='icrs')
            sep = coord.separation(match_coord).arcsec
            
            # Extract proper motion and parallax, handle masked values
            # Try uppercase and lowercase column names
            pmra = row.get('PMRA', row.get('pmra', np.nan))
            if np.ma.is_masked(pmra): pmra = np.nan
            
            pmdec = row.get('PMDEC', row.get('pmdec', np.nan))
            if np.ma.is_masked(pmdec): pmdec = np.nan
            
            # PLX might be PLX_VALUE or PLX
            plx = row.get('PLX_VALUE', row.get('PLX', row.get('plx_value', row.get('plx', np.nan))))
            if np.ma.is_masked(plx):
                plx = np.nan
            
            # MAIN_ID might be main_id
            main_id = row.get('MAIN_ID', row.get('main_id', 'Unknown'))
            
            # Parse Catalogue from ID (e.g., "Gaia DR2 ..." -> "Gaia")
            # Take the first word or known prefixes
            cat_name = 'Simbad'
            if isinstance(main_id, str):
                parts = main_id.split()
                if len(parts) > 0:
                    # Heuristics for common catalogues
                    if parts[0] in ['Gaia', 'Tycho', 'HIP', 'UCAC4', '2MASS', 'WISE', 'WISEP', 'AllWISE']:
                        cat_name = parts[0]
                    elif parts[0] == 'NAME':
                         # "NAME Barnard's Star" -> "Barnard's Star"? Or just keep full name?
                         # Maybe just use "Named Object" or the full name?
                         # Let's try to be smart. If it starts with NAME, use the rest.
                         if len(parts) > 1:
                             cat_name = " ".join(parts[1:])
                         else:
                             cat_name = main_id
                    else:
                        cat_name = parts[0]
            
            results.append({
                'ID': main_id,
                'RA': ra_val,
                'Dec': dec_val,
                'Distance_arcsec': sep,
                'Catalogue': cat_name,
                'PMRA': pmra,
                'PMDec': pmdec,
                'Plx': plx
            })
            
        return pd.DataFrame(results)
    except Exception as e:
        print(f"Error querying Simbad: {e}")
        return pd.DataFrame()

def plot_results(target_coord, results_df, radius_arcsec, output_file=None):
    """
    Plot the target and matched sources.
    """
    # Use WCS for proper axes
    # We create a simple WCS centered on the target
    # But for a simple scatter plot with small FOV, standard plotting with formatted ticks is easier and sufficient
    # as long as we format the ticks to HMS/DMS.
    
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Target
    ax.scatter(target_coord.ra.deg, target_coord.dec.deg, marker='+', s=300, c='black', label='Target', zorder=10)
    
    # Matches
    if not results_df.empty:
        # Define colors/markers for catalogues
        unique_cats = results_df['Catalogue'].unique()
        # Use a colormap that supports many categories or cycle
        colors = plt.cm.tab10(np.linspace(0, 1, len(unique_cats)))
        
        for i, cat in enumerate(unique_cats):
            group = results_df[results_df['Catalogue'] == cat]
            # Cycle markers if many catalogues?
            marker = 'o'
            
            ax.scatter(group['RA'], group['Dec'], label=cat, color=colors[i], marker=marker, s=50, alpha=0.8, zorder=5)
            
            for _, row in group.iterrows():
                # Draw faint line to target
                ax.plot([target_coord.ra.deg, row['RA']], [target_coord.dec.deg, row['Dec']], 
                        color=colors[i], linestyle='-', alpha=0.3, linewidth=0.5)
                
                # Annotate distance
                mid_ra = (target_coord.ra.deg + row['RA']) / 2
                mid_dec = (target_coord.dec.deg + row['Dec']) / 2
                ax.text(mid_ra, mid_dec, f"{row['Distance_arcsec']:.1f}\"", 
                        color=colors[i], fontsize=8, alpha=0.8, ha='center', va='center')
                
                # Proper Motion Correction
                if not pd.isna(row['PMRA']) and not pd.isna(row['PMDec']):
                    try:
                        # Propagate to now (approximate current time if not specified, let's use 2025.0 as per prop_motion default or current year)
                        # prop_motion function prints stuff, we might want to suppress it or just let it print.
                        # We need to create a SkyCoord for the source
                        source_coord = SkyCoord(ra=row['RA']*u.deg, dec=row['Dec']*u.deg, frame='icrs')
                        
                        # Distance for parallax
                        dist_pc = None
                        if not pd.isna(row['Plx']) and row['Plx'] > 0:
                            dist_pc = 1000.0 / row['Plx']
                        
                        # Current time
                        now_time = Time.now()
                        
                        # Capture stdout to suppress print from prop_motion if desired, or just let it run
                        # We'll just run it.
                        # Note: prop_motion expects strings for time usually, or Time objects? 
                        # Looking at prop_motion.py: old_time='J2000', new_time='2025-01-01...'
                        # It uses Time(old_time).
                        
                        # Simbad coords are usually J2000.
                        
                        # We need to temporarily redirect stdout to avoid cluttering the plot output if prop_motion is chatty
                        # But prop_motion is imported.
                        
                        # Let's use the logic from prop_motion directly or call it.
                        # Calling it is safer to reuse logic.
                        
                        # Suppress print
                        with open(os.devnull, 'w') as f, contextlib.redirect_stdout(f):
                             new_pos = prop_motion(source_coord, row['PMRA'], row['PMDec'], 
                                                  old_time='J2000', new_time=now_time.isot, 
                                                  distance=dist_pc, degrees=True)
                        
                        # Plot PM vector
                        ax.plot([row['RA'], new_pos.ra.deg], [row['Dec'], new_pos.dec.deg], 
                                color=colors[i], linestyle='--', alpha=0.6, linewidth=1)
                        ax.scatter(new_pos.ra.deg, new_pos.dec.deg, marker='x', color=colors[i], s=30, alpha=0.6)
                        
                    except Exception as e:
                        # print(f"PM Error: {e}")
                        pass

    # Search radius circle
    circle = plt.Circle((target_coord.ra.deg, target_coord.dec.deg), radius_arcsec/3600, color='gray', fill=False, linestyle='--', label=f'{radius_arcsec}" Radius')
    ax.add_patch(circle)
    
    # Formatting Axes to HMS/DMS
    # We can use a custom formatter
    
    def deg_to_hms(x, pos):
        h = int(x / 15)
        m = int((x / 15 - h) * 60)
        s = ((x / 15 - h) * 60 - m) * 60
        return f'{h:02d}h{m:02d}m{s:04.1f}s'

    def deg_to_dms(x, pos):
        d = int(x)
        m = int(abs(x - d) * 60)
        s = (abs(x - d) * 60 - m) * 60
        return f'{d:+03d}d{m:02d}m{s:04.1f}s'

    # Only apply if range is small enough to warrant it, but user asked for it.
    # However, standard matplotlib ticker might be messy with long labels.
    # Let's try to be smart. If the range is very small (arcseconds), we might just want offsets or careful formatting.
    # The user asked for "units of hms and dms".
    
    # We will use the formatters but maybe rotate labels.
    from matplotlib.ticker import FuncFormatter, MultipleLocator, AutoMinorLocator
    ax.xaxis.set_major_formatter(FuncFormatter(deg_to_hms))
    ax.yaxis.set_major_formatter(FuncFormatter(deg_to_dms))
    
    # Set ticks every 15 arcseconds (approx)
    # 15 arcsec = 30/3600 degrees

    tick_interval = 15.0 / 3600.0
    ax.xaxis.set_major_locator(MultipleLocator(tick_interval))
    ax.yaxis.set_major_locator(MultipleLocator(tick_interval))
    
    # Minor ticks: every 1 arcsecond?
    # 1 arcsec = 1/3600 degrees
    # If major is 30, minor could be 5 or 10?
    # AutoMinorLocator(5) puts 4 minor ticks (interval = major/5 = 6 arcsec).
    # Let's try AutoMinorLocator to keep it clean.
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    
    plt.xticks(rotation=45)
    
    plt.xlabel('RA (J2000)')
    plt.ylabel('Dec (J2000)')
    plt.title(f'Cone Search Results (Radius: {radius_arcsec}")')
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, which='major', linestyle=':', alpha=0.6)
    plt.grid(True, which='minor', linestyle=':', alpha=0.2)
    
    # Invert RA axis
    ax.invert_xaxis()
    
    # Adjust limits
    margin = radius_arcsec / 3600 * 1.2
    ax.set_xlim(target_coord.ra.deg + margin, target_coord.ra.deg - margin) # Inverted
    ax.set_ylim(target_coord.dec.deg - margin, target_coord.dec.deg + margin)
    
    plt.tight_layout()
    
    if output_file:
        plt.savefig(output_file)
        print(f"Plot saved to {output_file}")
    else:
        plt.show()

import contextlib # For suppressing stdout

@click.command()
@click.argument('coords', type=str)
@click.option('--radius', '-r', default=5.0, help='Search radius in arcseconds (default: 5).')
@click.option('--matches', '-m', default=5, help='Number of top matches to display (default: 5).')
@click.option('--radio-only', is_flag=True, default=False, help='Search only specific radio catalogues.')
@click.option('--save-csv', is_flag=True, default=False, help='Save results to a CSV file.')
@click.option('--plot/--no-plot', default=True, help='Generate a plot of the results.')
@click.option('--output-plot', default='search_plot.png', help='Filename for the plot if saved.')
def main(coords, radius, matches, radio_only, save_csv, plot, output_plot):
    """
    Perform a cone search for astronomical sources.
    
    COORDS should be a string with RA and Dec (e.g., "12:34:56 -45:67:89" or "188.73 -45.13").
    """
    
    try:
        # Suppress unicoord output
        with open(os.devnull, 'w') as f, contextlib.redirect_stdout(f):
             pos_eq, _ = unicoord(coords, galactic=False, display=True)
    except Exception as e:
        print(f"Error parsing coordinates: {e}")
        return

    print(f"Searching around {pos_eq.to_string('hmsdms')} with radius {radius}\"...")

    if radio_only:
        print("Performing Radio-only search...")
        df = perform_radio_search(pos_eq, radius)
    else:
        print("Performing Simbad search...")
        df = perform_simbad_search(pos_eq, radius)

    if df.empty:
        print("No matches found.")
        if plot:
             plot_results(pos_eq, df, radius, output_file=output_plot)
        return

    df['Distance_arcsec'] = df['Distance_arcsec'].astype(float)
    df = df.sort_values('Distance_arcsec').reset_index(drop=True)
    
    # Select columns to display
    disp_cols = ['ID', 'RA', 'Dec', 'Distance_arcsec', 'Catalogue']
    if 'PMRA' in df.columns:
        disp_cols.extend(['PMRA', 'PMDec'])
        
    print(f"\nTop {matches} Matches:")
    print(df[disp_cols].head(matches).to_string(index=False))
    
    if save_csv:
        filename = f"search_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
        df.to_csv(filename, index=False)
        print(f"\nResults saved to {filename}")

    if plot:
        plot_results(pos_eq, df, radius, output_file=output_plot)

if __name__ == '__main__':
    main()

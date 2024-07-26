from matplotlib import pyplot as plt
import numpy as np
import os
import itertools


def main():
    directory = 'C:/Users/khoaw/OneDrive/Personal Documents/UT - REU/Data/QE Data/07222024' # TO CHANGE

    # Initialize program variables
    meas_name = {
        'cell_type': '',
        'sample_num': '',
        'cell_num': '',
        'temp_num': 0
    }
    temp_list = []
    Jsc_AM15G_list = []
    Jsc_AM0_list = []
    Eg_list = []
    legend_list = []
    last_file = ''

    # Store list of different plot styles
    color = itertools.cycle(('b', 'g', 'r', 'c', 'm', 'y', 'k'))
    line_style = itertools.cycle(('-', ':', '--', '-.'))

    for i, files in enumerate(os.listdir(directory)):
        if not files.startswith('QE') and files == 'Plots':
            continue

        # If new sample, plot parameters vs temp and clear all lists
        file_name = files.split('-')
        if i != 1 and file_name[0] != last_file:
            plot_param_vs_temp(temp_list, Jsc_AM15G_list, Jsc_AM0_list, Eg_list, meas_name, True)
            plt.savefig(directory+'/Plots/' + last_file + '_trends_best.png', bbox_inches='tight')
            plot_param_vs_temp(temp_list, Jsc_AM15G_list, Jsc_AM0_list, Eg_list, meas_name, False)
            plt.savefig(directory + '/Plots/' + last_file + '_trends_conn.png', bbox_inches='tight')

            plt.clf()
            temp_list.clear()
            Jsc_AM15G_list.clear()
            Jsc_AM0_list.clear()
            Eg_list.clear()
            legend_list.clear()

        # Open QE Measurement files
        wavelength_list, QE_list = np.loadtxt(os.path.join(directory, files), delimiter='\t',
                                              skiprows=17, unpack=True, usecols=(0, 1))

        # Get Jsc parameters
        Jsc_AM15G, Jsc_AM0, Eg = run_calculations(wavelength_list, QE_list)

        # Get sample ID from file name (MAY HAVE TO CHANGE)
        sample_IDs = file_name[0].split('_')
        meas_name['cell_type'] = sample_IDs[1]
        meas_name['sample_num'] = sample_IDs[2]
        meas_name['cell_num'] = sample_IDs[3]
        meas_name['temp_num'] = int(file_name[1][:len(file_name[1])-1])

        # Add to lists
        temp_list.append(meas_name['temp_num'])
        Jsc_AM15G_list.append(Jsc_AM15G)
        Jsc_AM0_list.append(Jsc_AM0)
        Eg_list.append(Eg)

        # Plot the data
        plt.plot(wavelength_list, QE_list, color=next(color), linestyle=next(line_style))

        # Add labels to the plot
        plt.title('QE of {} PVSK - Sample {} Cell {}'.format(
            meas_name['cell_type'], meas_name['sample_num'], meas_name['cell_num']))
        plt.xlabel('Wavelength [nm]')
        plt.ylabel('QE [%]')
        plt.ylim(0, 100)

        # Add to the legend
        legend_list.append(str(meas_name['temp_num']) + 'C')
        plt.legend(legend_list, loc='upper right')

        # Save the plot
        if not os.path.exists(directory + '/Plots'):
            os.makedirs(directory + '/Plots')
        plt.savefig(directory+'/Plots/'+file_name[0] + '.png', bbox_inches='tight')

        # Remember the last file name
        last_file = file_name[0]

    # Plot last sample
    plot_param_vs_temp(temp_list, Jsc_AM15G_list, Jsc_AM0_list, Eg_list, meas_name, True)
    plt.savefig(directory + '/Plots/' + last_file + '_trends_best.png', bbox_inches='tight')
    plot_param_vs_temp(temp_list, Jsc_AM15G_list, Jsc_AM0_list, Eg_list, meas_name, False)
    plt.savefig(directory + '/Plots/' + last_file + '_trends_conn.png', bbox_inches='tight')


def run_calculations(wavelength_list, QE_list):
    # Constants
    h = 6.626e-34  # Planck's constant
    q = 1.602e-19  # Charge of an electron
    c = 3e8  # Speed of light

    # Get AM1.5G and AM0 irradiance values with corresponding wavelengths
    wavelength_am15, irradiance_am15 = np.loadtxt('am1.5g-PVLighthouse.txt', unpack=True)
    wavelength_am0, irradiance_am0 = np.loadtxt('am0-PVLighthouse.txt', unpack=True)

    # Create lists to store the Jsc operands data
    integrand_AM15G = []
    integrand_AM0 = []
    Jsc_am15g_sum = 0
    Jsc_am0_sum = 0
    QE_der_list = []
    last_QE = 0

    for i, (wave_len, QE_wave_len) in enumerate(zip(wavelength_list, QE_list)):
        # Method 1 Calculations
        AM15G_irrad = irradiance_am15[np.where(wavelength_am15 == wave_len)[0][0]]  # Find corresponding H
        AM0_irrad = irradiance_am0[np.where(wavelength_am0 == wave_len)[0][0]]
        integrand_AM15G.append((QE_wave_len / 100) * (AM15G_irrad * wave_len))  # EQE(wavelen) * N(wavelen)
        integrand_AM0.append((QE_wave_len / 100) * (AM0_irrad * wave_len))

        # Method 2 Calculations
        if i == 0:
            Jsc_am15g_sum = Jsc_am15g_sum + (QE_wave_len / 100) * ((AM15G_irrad * wave_len) * 1e-13) * q / (h * c)
            Jsc_am0_sum = Jsc_am0_sum + (QE_wave_len / 100) * ((AM15G_irrad * wave_len) * 1e-13) * q / (h * c)
        else:
            Jsc_am15g_sum = Jsc_am15g_sum + (QE_wave_len / 100) * ((AM15G_irrad * wave_len) * 1e-13) * 10 * q / (h * c)
            Jsc_am0_sum = Jsc_am0_sum + (QE_wave_len / 100) * ((AM15G_irrad * wave_len) * 1e-13) * 10 * q / (h * c)

        # Find the QE differentials for Eg calculation
        if i != 0:
            QE_der_list.append(last_QE - QE_wave_len)
        last_QE = QE_wave_len

    # Integrate over wavelengths in meters and convert to mA/cm^2
    wave_met = [i * 1e-9 for i in wavelength_list]
    Jsc_AM15G = np.trapezoid(integrand_AM15G, wave_met) * q / (h * c) / 10
    Jsc_AM0 = np.trapezoid(integrand_AM0, wave_met) * q / (h * c) / 10

    # Calculate band gap
    Eg = 1240 / (wavelength_list[QE_der_list.index(max(QE_der_list))]+5.)

    return Jsc_AM15G, Jsc_AM0, Eg


def plot_param_vs_temp(temp_list, Jsc_AM15G_list, Jsc_AM0_list, Eg_list, meas_name, flag):
    plt.clf()
    fig, ax = plt.subplots()
    fig.subplots_adjust(right=0.75)

    # Create new y-axes
    twin1 = ax.twinx()
    twin2 = ax.twinx()

    # Move extra y-axis spines to the right
    twin2.spines.right.set_position(("axes", 1.2))

    # Plot points with temperature coefficients in the label
    p1, = ax.plot(temp_list, Jsc_AM15G_list, label="Jsc_AM1.5G ({:.4f} mA/cm{}/C)".format(
        np.polyfit(temp_list, Jsc_AM15G_list, 1)[0], get_super('2')),
        marker='+', linestyle='' if flag else '-', color='r')
    p2, = twin1.plot(temp_list, Jsc_AM0_list, label="Jsc_AM0 ({:.4f} mA/cm{}/C)".format(
        np.polyfit(temp_list, Jsc_AM0_list, 1)[0], get_super('2')),
        marker='.', linestyle='' if flag else '-', color='b')
    p3, = twin2.plot(temp_list, Eg_list, label="Band-Gap ({:.4f} eV/C)".format(
        np.polyfit(temp_list, Eg_list, 1)[0]),
        marker='*', linestyle='' if flag else '-', color='y')

    # Plot lines of best fit
    if flag:
        ax.plot(np.unique(temp_list),
                np.poly1d(np.polyfit(temp_list, Jsc_AM15G_list, 1))(np.unique(temp_list)), color='r')
        twin1.plot(np.unique(temp_list),
                   np.poly1d(np.polyfit(temp_list, Jsc_AM0_list, 1))(np.unique(temp_list)), color='b')
        twin2.plot(np.unique(temp_list),
                   np.poly1d(np.polyfit(temp_list, Eg_list, 1))(np.unique(temp_list)), color='y')

    # Add titles to axes
    plt.title('Jsc and Eg of {} PVSK - Sample {} Cell {} at Varying Temperatures'.format(
        meas_name['cell_type'], meas_name['sample_num'], meas_name['cell_num']))
    ax.set(xlim=(14, 46), xlabel="Temperature [C]", ylabel="Jsc AM1.5G [mA/cm{}]".format(get_super('2')))
    twin1.set(ylabel="Jsc_AM0 [mA/cm{}]".format(get_super('2')))
    twin2.set(ylabel="Band-Gap [eV]")

    # Set color of axis titles
    ax.yaxis.label.set_color(p1.get_color())
    twin1.yaxis.label.set_color(p2.get_color())
    twin2.yaxis.label.set_color(p3.get_color())

    # Set color of axis labels
    ax.tick_params(axis='y', colors=p1.get_color())
    twin1.tick_params(axis='y', colors=p2.get_color())
    twin2.tick_params(axis='y', colors=p3.get_color())

    # Add legend to graph of plot labels
    ax.legend(handles=[p1, p2, p3])


def get_super(x):
    normal = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+-=()"
    super_s = "ᴬᴮᶜᴰᴱᶠᴳᴴᴵᴶᴷᴸᴹᴺᴼᴾQᴿˢᵀᵁⱽᵂˣʸᶻᵃᵇᶜᵈᵉᶠᵍʰᶦʲᵏˡᵐⁿᵒᵖ۹ʳˢᵗᵘᵛʷˣʸᶻ⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻⁼⁽⁾"
    res = x.maketrans(''.join(normal), ''.join(super_s))
    return x.translate(res)


if __name__ == "__main__":
    main()

'''
This is a library to support processing, analising and plotting the data
form the single-shot qubit readout experiment (with a postselection).


Vladimir Milchakov
vladimir_ph@protonmail.com
'''



import numpy as np
from matplotlib import pyplot as plt

################################################################################
### My Design ###
my_colors_dict = {  'redberry'      :'#970000',
                    'tamarillo'     :'#8b1212',
                    'venetianred'   :'#770023', ## main
                    'monarch'       :'#850731',
                    'toledo'        :'#40001b',
                    'shipgrey'      :'#423B4B',
                    'charde'        :'#20202c',
                    'woodsmoke'     :'#171719',
                    'mediumpurple'  :'#795FD7',
                    'curciousblue'  :'#2f99e2',
                    'electric'      :'#795fd7',
                    'deus_ex_gold'  :'#eda723',
                    'meduza'        :'#b5895b',
                    'meduza_gold'   :'#B5965B',
                    'blob_g'        :'#546ea9',
                    'blob_e'        :'#fd6e6e',
                    'blob_post'     :'#1b346e',
                    'g_state_mark'  :'#849cd4',
                    # 'g_state_mark'  :'#2b488a',
                    'e_state_mark'  :'#fe9393',
                    # 'e_state_mark'  :'#a53030',
                    'meduza_dark'   :'#262626',
                    'gauss_green'   :'#bacc5f'

}

dark_scheme = {     'blob_g'            :my_colors_dict['blob_g'],
                    'blob_e'            :my_colors_dict['blob_e'],
                    'blob_prepulse'     :'g',
                    'color_g_mean'      :my_colors_dict['g_state_mark'],
                    'color_e_mean'      :my_colors_dict['e_state_mark'],
                    'color_g_cross'     :my_colors_dict['deus_ex_gold'],
                    'color_e_cross'     :my_colors_dict['deus_ex_gold'],
                    'color_dist'        :my_colors_dict['deus_ex_gold'],
                    'color_g_vector'    :my_colors_dict['g_state_mark'],
                    'color_e_vector'    :my_colors_dict['e_state_mark'],
                    'color_void'        :my_colors_dict['deus_ex_gold'],
                    'color_zero'        :my_colors_dict['meduza_gold'],
                    'color_zero_vector' :my_colors_dict['meduza_gold'],
                    'fig_face_color'    :my_colors_dict['meduza_dark'] ,
                    'fig_border_color'  :'r',
                    'bg_color'          :'k',
                    'grid_color'        :my_colors_dict['meduza_gold'],
                    'title_color'       :my_colors_dict['meduza_gold'],
                    'legend_color'      :my_colors_dict['meduza_dark'],
                    'legend_text_color' :my_colors_dict['meduza_gold'],
                    'legend_frame_color':my_colors_dict['meduza_gold'],
                    'AXES_COLOR'        :my_colors_dict['meduza_gold'],
                    'legend_alpha'       :0.7,
                    'grid_transp'        :0.5,
                    'lw_cross'           :0.8,
                    'vector_bw_blobs'    :1,
                    'vector_state_lw'    :1,
                    'vector_zero_lw'     :1,
                    'point_scatter_size' :1,
                    'onplot_mark_size'   :10


}

bright_scheme = {   'blob_g'            :my_colors_dict['blob_g'],
                    'blob_e'            :my_colors_dict['blob_e'],
                    'blob_prepulse'     :'g',
                    'color_g_mean'      :my_colors_dict['g_state_mark'],
                    'color_e_mean'      :my_colors_dict['e_state_mark'],
                    'color_g_cross'     :'k',
                    'color_e_cross'     :'k',
                    'color_dist'        :'k',
                    'color_g_vector'    :'b',
                    'color_e_vector'    :'r',
                    'color_void'        :'lightgrey',
                    'color_zero'        :'lightgrey',
                    'color_zero_vector' :'lightgrey',
                    'fig_face_color'    :'white' ,
                    'fig_border_color'  :'r',
                    'bg_color'          :'white',
                    'grid_color'        :'lightgrey',
                    'title_color'       :'k',
                    'legend_color'      :'lightgrey',
                    'legend_text_color' :'k',
                    'legend_frame_color':'lightgrey',
                    'AXES_COLOR'        :'k',
                    'legend_alpha'       :0.8,
                    'grid_transp'        :0.7,
                    'lw_cross'           :1,
                    'vector_bw_blobs'    :1,
                    'vector_state_lw'    :1,
                    'vector_zero_lw'     :1

}




def set_font(filename='Forza-Book.ttf'):
    '''
    Takes the name of file.ttf
    File must be in folder of
    Python27\Lib\site-packages\matplotlib\mpl-data\fonts\ttf
    return prop
    wich could be use as in an example:
    ax.set_title('This is a special font: {}'.format(fname), fontproperties=prop)
    fontproperties = set_font('Forza-Book.ttf')
    '''
    ## special font from ttf
    ###+===============================######
    import os
    from matplotlib import font_manager as fm, rcParams

    str_adress_name = "fonts/ttf/" + filename

    fpath = os.path.join(rcParams["datapath"], str_adress_name)
    prop = fm.FontProperties(fname=fpath)
    # fname = os.path.split(fpath)[1]

    return prop

################################################################################

################################################################################
### This needs only for fit the hists (maybe future) ###
def sameside(a,b, ref=0):
    '''
    if ref = 0
    return True if a and b has a same sign
    for different ref return True if a and b on the one side from ref
    '''
    if (a == ref) or (b==ref):
        return True
    a_sign = a > ref
    b_sign = b > ref
    if a_sign == b_sign:
        return True
    else:
        return False

def where_sameside(a_arr, b, ref=0):
    '''
    a_arr = array, b =value_to_compare_with, and ref = reference value
    Returns array of booleans:
    True - on positions where value of both arrays have same side to ref
    False where side is different
    '''
    result_array = np.zeros(len(a_arr))
    for i in range(len(a_arr)):
        result_array[i] = sameside(a_arr[i] ,b, ref=ref)

    return result_array

################################################################################

################################################################################
### math functions ###
def my_stround(num, digits=3, withminus=False,toright=False):
    '''
    Takes float. Return string. With given number of digits.
    '''
    digits = abs(int(digits))
    dig_remain = digits

    ### work with first char
    if num >= 0:
        if withminus:
            minus = ' '
            dig_remain = dig_remain -1
        else:
            minus = ''
    else:
        minus = '-'
        dig_remain = dig_remain -1

    num = abs(num)

    ### separate two parts
    left_part = int(num)
    right_part = num - left_part

    dig_remain = dig_remain - len(str(left_part))

    if dig_remain < 1:
        string = minus + str(left_part)
    elif dig_remain ==1:
        if toright:
            string = ' ' + minus + str(left_part)
        else:
             string = minus + str(left_part) + ' '
    else:
        right_part = round(right_part, dig_remain-1)
        string = minus + str( left_part + right_part  )

    if len(string) < digits:
        for i in range(  digits-len(string)  ):
            string = string + '0'

    return string

def complex_num_relationships(re1,im1,re2,im2):
    '''
    return distance and angle between two complex numbers
    tekes re and im of both numbers
    '''
    c1 = re1 + 1j*im1
    c2 = re2 + 1j*im2
    theta = np.angle(c2-c1, deg=True)
    distance = np.sqrt( (re2-re1)**2 + (im2-im1)**2  )

    return [distance, theta]

def angle_three_points(x1,y1,x2,y2,x3,y3):
    '''
    Find an angle between three points.
    Takes 6 coordinats <ABC
    return angle in grad
    '''

    [dist, theta12] = complex_num_relationships(x2,y2,x1,y1)
    if theta12 < 0:
        theta12 = 360 + theta12
    [dist, theta23] = complex_num_relationships(x3,y3,x2,y2)
    if theta23 < 0:
        theta23 = 360 + theta23

    angle = theta23 - theta12 - 180

     # to be always in range [-360:360]
    if angle> 0:
        angle = angle % 360
    else:
        angle = -angle % 360
        angle = -angle

    return angle

def change_basis_point(x,y, x0, y0, theta):
    '''
    Change the basis of one point.
    takes coordinates of one point and bassis as [x0,y0, theta(grad)]
    return coordinatees in a new basis
    First shift, than rotate. Save shift as if shift was after rotation
    '''
    ### shift the point
    x1 = x -x0
    y1 = y -y0
    #convert theta to radians
    theta = np.radians(theta)
    ### rotate the vector clockwise
    x2 =  x1*np.cos(theta) + y1*np.sin(theta)
    y2 = -x1*np.sin(theta) + y1*np.cos(theta)
    ########

    return [x2,y2]

def change_basis_blobs_inf(x0, y0, theta, *args_datapairs):
    '''
    change basis of all sequence
    takes x0,y0 and theta of a new basis
    AND  few of ndarrays of points in shape [x0_arr, y0_arr], [x1_arr, y1_arr], ...[xn_arr, yn_arr]
    returns few of ndarrays of points in New Basis in shape of list of lists of arrays: [  [x0_arr, y0_arr], [...], ... [xn_arr, yn_arr] ]
    (theta in degrees)
    '''
    result_list = []

    for data_pair in args_datapairs:
        data_x = data_pair[0]
        data_y = data_pair[1]
        if len(data_x) != len(data_y):
            print '___WARNING change_basis_blobs_inf() not the same length of arrays!'

        #for each given sequence - do the basis changing
        data_x_norm = np.ones_like(data_x)  #prepare empty array
        data_y_norm = np.ones_like(data_y)

        # for i in range(len(re_g)):
        for i in range(len(data_x)):
            [ data_x_norm[i], data_y_norm[i] ] = change_basis_point(data_x[i], data_y[i], x0, y0, theta)    #g state

        data_norm = [data_x_norm, data_y_norm]
        result_list.append(data_norm)   #add sequence in a new basis to the list

    return result_list

def centers_two_blobs(re_g, im_g, re_e, im_e):
    '''
    Searching the centers of two blobs on i-q plane
    Takes ndarrays of Re_g, Im_g, Re_e, Im_e
    return coordinates of center list = [c_re_g, c_im_g, c_re_e, c_im_e]
    '''
    c_re_g = np.mean(re_g)
    c_im_g = np.mean(im_g)
    c_re_e = np.mean(re_e)
    c_im_e = np.mean(im_e)

    return [c_re_g, c_im_g, c_re_e, c_im_e ]

def crop_fluctuations(re_g, im_g, re_e,im_e, void_re, void_im, coeff_shift = 2.0, crop=True):
    '''
    takes re_g, im_g, re_e, im_e data
    finds the main area
    returns values of limits [leftlim, rightlim, toplim, bottomlim]
    coeff_shift - how many std-values shift from mean-value
    '''
    if crop:
        ## for X axis
        re_g_m = np.mean(re_g)     #mean value of g
        re_g_d = np.std(re_g)       #dispersion value of g
        re_e_m = np.mean(re_e)
        re_e_d = np.std(re_e)
        ## for Y axis
        im_g_m = np.mean(im_g)     #mean value of g
        im_g_d = np.std(im_g)       #dispersion value of g
        im_e_m = np.mean(im_e)
        im_e_d = np.std(im_e)
        ### possible X(re) limits
        re1 = re_g_m + re_g_d * coeff_shift
        re2 = re_g_m - re_g_d * coeff_shift
        re3 = re_e_m + re_e_d * coeff_shift
        re4 = re_e_m - re_e_d * coeff_shift
        ### possible Y(im) limits
        im1 = im_g_m + im_g_d * coeff_shift
        im2 = im_g_m - im_g_d * coeff_shift
        im3 = im_e_m + im_e_d * coeff_shift
        im4 = im_e_m - im_e_d * coeff_shift

        leftlim = np.min([ re1,re2,re3,re4, 0, void_re])
        rightlim = np.max([ re1,re2,re3,re4, 0, void_re ])
        toplim = np.max([ im1,im2,im3,im4, 0, void_im ])
        bottomlim = np.min([ im1,im2,im3,im4, 0, void_im ])
    else:
        leftlim = np.min([ re_e, re_g, 0 ])
        rightlim = np.max([ re_e, re_g, 0 ])
        toplim = np.max([ im_e, im_g, 0 ])
        bottomlim = np.min([ im_e, im_g, 0 ])

    return [leftlim, rightlim, toplim, bottomlim]

def get_count_states(x_g, x_e, threshold):
    '''
    Function to count number of values according to threshold
    Returns dictionary
    n_left_g = number elements x_g less than threshold
    n_left_e = number elements x_e less than threshold
    n_right_g = number elements x_g more than threshold
    n_right_e = number elements x_e more than threshold
    p_left_g -probabilities
    p_left_e
    p_right_g
    p_right_e
    sighn_ge - +1 if g less than e; or -1 if g more than e
    p_ge = P to measure e when g was prepared
    '''
    def e_morethan_g(p_left_g, p_left_e, p_right_g, p_right_e):
        '''
        Function for define where is g state and where is e state relative to the threshold
        Takes probabilities of g and e states to be on the left and right parts
        Returns True if |e> > |g>, and False if |g> > |e>
        '''
        if (p_left_g > p_right_g) and (p_left_e < p_right_e):
            return True
        elif (p_left_g < p_right_g) and (p_left_e > p_right_e):
            return False
        elif (p_left_g > p_right_g) and (p_left_e > p_right_e):
            if p_left_g > p_left_e:
                return True
            else:
                return False
        elif (p_left_g < p_right_g) and (p_left_e < p_right_e):
            if p_left_g > p_left_e:
                return True
            else:
                return False

    if ( len(x_g) <1 ) or ( len(x_e) <1 ):
            print 'error get_count_states. No data'
            return None

    #calculate how many g points are less or more than threshold value
    n_left_g  = 0
    n_right_g = 0
    for val in x_g:
        if val < threshold:
            n_left_g  +=1
        else:
            n_right_g +=1
    #calculate how many e points are less or more than threshold value
    n_left_e  = 0
    n_right_e = 0
    for val in x_e:
        if val < threshold:
            n_left_e  +=1
        else:
            n_right_e +=1
    #probabilities of g or e states been bigger or smaller than threshold
    #### DO YOU REALLY NEED PROBABILITIES NOW?

    p_left_g = 1.0*n_left_g / len(x_g)
    p_left_e = 1.0*n_left_e / len(x_e)
    p_right_g = 1.0*n_right_g / len(x_g)
    p_right_e = 1.0*n_right_e / len(x_e)


    ### sign ###
    # check where is e and where is g. Define P_ij
    if e_morethan_g(p_left_g, p_left_e, p_right_g, p_right_e):
        sign = +1
        p_gg = p_left_g
        p_ee = p_right_e
        p_ge = p_left_e # P to measure e when g was prepared
        p_eg = p_right_g # P to measure g when e was prepared
    else:
        sign = -1
        p_gg = p_left_e
        p_ee = p_right_g
        p_ge = p_right_e # P to measure e when g was prepared
        p_eg = p_left_g # P to measure g when e was prepared

    dict_count = {
    'n_left_g'  :   n_left_g,
    'n_left_e'  :   n_left_e,
    'n_right_g' :   n_right_g,
    'n_right_e' :   n_right_e,
    'p_left_g'  :   p_left_g,
    'p_left_e'  :   p_left_e,
    'p_right_g' :   p_right_g,
    'p_right_e' :   p_right_e,
    'sign_ge'   :   sign,
    'p_gg'      :   p_gg,
    'p_ee'      :   p_ee,
    'p_ge'      :   p_ge,   # P to measure e when g was prepared
    'p_eg'      :   p_eg,   # P to measure g when e was prepared
    }

    return dict_count

def get_overlap_error(gauss_p_g, gauss_p_e, threshold, sign_ge):
    '''
    Function to calculate gaussian overlap errors and fidelity
    returns dictionary with errors_overlap and fidelities
    Arguments: gauss_p_g, gauss_p_e, threshold, sign_ge
    gauss_p_g  - parameters of gaussian fit of |g> - state
    gauss_p_e  - parameters of gaussian fit of |e> - state
    threshold  - threshold value
    sign_ge    - if |g> less than |e> it is [+1]. (otherwise [-1] 0
    '''
    import scipy

    def errfunc(x, x_c=0, sig=1.0):
        '''
        modification for standard error-function from scipy
        to integrate gaussian form -infinity to x
        And also to take into account sigma and x_center
        '''
        y = scipy.special.erf( (x - x_c)/sig )
        y = 0.5 + 0.5*y
        return y

    ##################################

    [x_center_g, sigma_g] = [ gauss_p_g[0], gauss_p_g[1] ]
    [x_center_e, sigma_e] = [ gauss_p_e[0], gauss_p_e[1] ]

    if sign_ge == +1:
        error_g = 1.0 -errfunc(threshold, x_c=x_center_g, sig=sigma_g)
        error_e =      errfunc(threshold, x_c=x_center_e, sig=sigma_e)

    elif sign_ge == -1:
        error_g =      errfunc(threshold, x_c=x_center_g, sig=sigma_g)
        error_e = 1.0 -errfunc(threshold, x_c=x_center_e, sig=sigma_e)

    else:
        print 'Error. Something wrong with sign_ge. It must be -1 or +1'
        error_g = 0
        error_e = 0

    f_g_over = 1.0 - error_g
    f_e_over = 1.0 - error_e
    f_over = 1.0 - 0.5*( error_g+ error_e )

    dict_overlap = {
        'err_o_g': error_g,
        'err_o_e': error_e,
        'f_o_g'  : f_g_over,
        'f_o_e'  : f_e_over,
        'f_o'  : f_over
    }

    return dict_overlap

def max_possible_fid(t_read_time, T1 = 3500.):
    '''
    Calculate maximum possible fidelity for given t_read and T1
    '''
    p_ee = np.exp( -t_read_time/T1)
    p_gg = np.exp( -t_read_time/T1)
    p_eg = 1.0 - p_ee
    p_ge = 1.0 - p_gg
    F = 1.0 - 0.5*( p_ge +p_eg )
    return F
################################################################################

################################################################################
class Histogram:
    '''
    A class to containt histograms and it's fit data

    Variables:
    hist        - histogram itself = tuple: [val_arr, bins_arr]
    hist_xy     - histogram ramake to len(val) = len(bins) [val_arr, x_arr]
    gauss_fit   - also tuple with [val_arr, x_arr] but for fitted curve
    gauss_param - parameters of fit []  #[ center_v, sigma_v, std_v, max_y, min_y ]

    Method:
    fit( threshold, crop_sign )
        crop_sign - +1 if you need to take part more than threshold
                    -1 otherwise
    '''
    hist        = None
    hist_xy     = None
    gauss_fit   = None
    gauss_param = None  #[ center_v, sigma_v, std_v, max_y, min_y ]

    #########################################

    def __init__(self, x, nbins=100):
        '''
        takes 1D array of data
        '''
        self.hist = np.histogram(x, bins=nbins)

        def hist_tuple_to_xy_func(h_tuple):
            '''
            plt.hist automaticly return a tuple (values_array, coordinate_array)
            and they have different lenghtes bacouse of this:   .-.-.-.   :(3 lines, 4 points)
            this function reshape the coordinate array to make possible to match values to cordinates
            ### Thin it out!
            '''
            if h_tuple is None:
                return None
            vals = h_tuple[0]
            cords = h_tuple[1]

            cords_1 = np.zeros_like(cords[:-1])
            for i in range(len(cords_1)):
                cords_1[i] = np.mean([  cords[i], cords[i+1] ])

            return [vals, cords_1]

        self.hist_xy = hist_tuple_to_xy_func(self.hist)

    def fit(self, threshold, crop_sign, do_crop=True):
        '''
        decorator for function fit_gauss
        fit itself
        '''
        def fit_gauss(hist_x_xy, crop_thr=None, crop_sign = +1):
            '''
            Takes histogram-tuple [x, vals]
            fits
            returns histogram-tuple [x, vals_fit] and list of gauss-parameters
            crop_thr: if it is None - than we fit all the data. If it is threshold_value than we do crop
            if crop_sign is positive - than we take points more than threshold
            if crop_sign is negative - than we take points less than threshold
            '''

            def get_only_one_state_histxy(histxy, th=0, sign_ge=+1, state=+1):
                '''
                This function return only the points which is above (below) threshold
                th is threshold value
                sign_ge: {+1 if e>g, -1 otherwise}
                state is:
                for g-state state = -1
                for e-state state = +1
                '''
                x = histxy[1]
                y = histxy[0]

                if sign_ge * state > 0:
                    ind = np.where(x < th)
                elif sign_ge * state < 0:
                    ind = np.where(x > th)
                else:
                    print 'error of input cut_one_state_from_x(). sign_ge and state can be +-1 only'

                x_crop = np.delete(x, ind)
                y_crop = np.delete(y, ind)
                hist_crop = [y_crop, x_crop]

                return hist_crop

            def hwhh_xy(x,y):
                '''
                This function takes x and y sequence
                returns half width at half height
                Works for gaussian type of function (or lorenzian etc.)
                '''
                half_heigh = np.max(y)/2.0
                center_x = x[np.argmax(y)]

                list_x_in = []
                for i in range(len(y)):
                    if y[i] > half_heigh:
                        list_x_in.append(x[i])

                x_max_edge = np.max(list_x_in)
                x_min_edge = np.min(list_x_in)

                hwhh = (x_max_edge - x_min_edge)/2.0

                return hwhh

            #################################################

            ### to crop or not to crop
            if crop_thr is None:
                y = hist_x_xy[0]
                x = hist_x_xy[1]
                x_axis = x
            else:
                if crop_sign == +1:
                    hist_x_xy_crop = get_only_one_state_histxy( hist_x_xy, th=crop_thr, sign_ge=+1, state=crop_sign )
                elif crop_sign == -1:
                    hist_x_xy_crop = get_only_one_state_histxy( hist_x_xy, th=crop_thr, sign_ge=+1, state=crop_sign )
                else:
                    print 'warning in fit_gauss(). crop_sign must be equal to +-1'
                    return None
                y = hist_x_xy_crop[0]
                x = hist_x_xy_crop[1]
                x_axis = hist_x_xy[1]


            ### make an expectatin
            expec_min_y    = 0
            expec_max_y    = np.max(y)
            expec_center_x = x[np.argmax(y)]
            expec_std_x    = hwhh_xy(x,y)
            gauss_p0 = [expec_min_y, expec_max_y, expec_center_x, expec_std_x]

            ### fit g state by single gaus
            import fit
            gauss = fit.Gaussian()
            gauss.set_data(x, y)
            gauss_p = gauss.fit(gauss_p0, fixed=[0])
            y_gausfit = gauss.func(gauss_p)

            # remake the gaussian again with given parameters, but for all x-axis
            if crop_thr is not None:
                gauss1 = fit.Gaussian()
                gauss1.set_data(x_axis, y)
                y_gausfit = gauss1.func(gauss_p)


            ### here we already know exact centers!
            min_y       = gauss_p[0]
            max_y       = gauss_p[1]
            center_v    = gauss_p[2]
            fw_at_06    = gauss_p[3]

            sigma_v  = fw_at_06/2

            h_x_fit     = [ y_gausfit, x_axis ]
            p_list      = [ center_v, sigma_v, fw_at_06, max_y, min_y ]
            return [h_x_fit, p_list]

        if do_crop:
            [ hist_x_fit, gaus_par_x ]  = fit_gauss( self.hist_xy, crop_thr=threshold, crop_sign=crop_sign )
        else:
            [ hist_x_fit, gaus_par_x ]  = fit_gauss( self.hist_xy)

        self.gauss_fit           = hist_x_fit
        self.gauss_param         = gaus_par_x

        return True


class SSResult:
    '''
    This class is linked with one single-shot measurements.
    Data takes only from file during inizialization.
    Contains re and im parts of measured g and e states and also of two pre-measurements
    '''
    ###---------------------------###
    ###### DATA ######
    ### Raw Data ###
    CONVERT_TOMV=True       #Flag of convertion from [V] to [mV]

    void_re = 0
    void_im = 0
    re_g     = np.array([])
    im_g     = np.array([])
    re_e     = np.array([])
    im_e     = np.array([])
    re_g_pre = np.array([])
    im_g_pre = np.array([])
    re_e_pre = np.array([])
    im_e_pre = np.array([])

    ### Normalised Data ###
    void_x = 0
    void_y = 0
    x_g         = None
    x_e         = None
    x_g_pre     = None
    x_e_pre     = None
    y_g         = None
    y_e         = None
    y_g_p       = None
    y_e_pre     = None
    ### Save parameters of transformation of data
    THETA = 0       ## angle of rotation to go from re-im to x-y
    SHIFT_X = 0     ## x shift was done after rotation
    SHIFT_Y = 0     ## y shift was done after rotation

    ### Normalised Data after Postselection ###
    x_g_select  = None
    x_e_select  = None
    y_g_select  = None
    y_e_select  = None

    ## Histograms ### (it is all objects of class Histograms)
    hist_x_g         = None
    hist_x_e         = None
    hist_x_g_pre     = None
    hist_x_e_pre     = None
    hist_x_g_select  = None
    hist_x_e_select  = None

    hist_y_g         = None
    hist_y_e         = None

    ###---------------------------###


    ######################################################
    ####  MetaData -----------------##
    ######  ######

    ### Time ###
    datafile = ''
    timestamp = ''

    ### Parameters ###
    paramfile = ''
    # dict_param = {
    # 'freq_read':0, 'power1':0, 't_read':0, 'rudat':0, 'freq_q':0,
    #  'power2':0, 'rudat2':0, 'tpi':0, 'nsigma':0, 'cur':0
    #  }
    dict_param = None
    ################-----------------##
    ######################################################


    ######################################################
    ##### CENTERS, THRESHOLD -----------------##
    ### e-g-State definition ###
    threshold = None
    sign_ge = +1     #COULD BE +-1 depends on e>g or g>e

    ### centers of blobs (in first approximation - np.mean, after taken from gauss fit)
    center_x_g = None
    center_x_e = None
    center_x_g_select = None # it is different. I dont know why yet
    center_x_e_select = None
    sizeblob_x_g = None
    sizeblob_x_e = None
    sizeblob_y_g = None
    sizeblob_y_e = None
    squizingblob_g = None
    squizingblob_e = None
    ### crosses of size of each blolbs
    ## frame_g_reim = [ax,ay,bx,by,cx,cy,dx,dy]
    frame_g_xy = None
    frame_e_xy = None

    # ### Back_rotated real centers of blobs, extracted from histograms fits
    center_re_g = None
    center_re_e = None
    center_im_g = None
    center_im_e = None
    ### crosses of size of each blolbs in real raw axes
    ## frame_g_reim = [ax,ay,bx,by,cx,cy,dx,dy]
    frame_g_reim = None
    frame_e_reim = None

    ################-----------------##
    ######################################################

    ######################################################
    ### dimensionless R-axis (like X in [V])#####-------##
    r_g         = None
    r_e         = None
    r_g_pre     = None
    r_e_pre     = None
    r_g_select  = None
    r_e_select  = None

    center_r_g  = None
    center_r_e  = None
    void_r      = None
    void_r      = None
    threshold_r = None
    center_r_g_select = None
    center_r_e_select = None

    ### parameter S = (2/sigma_hist)^2 - 'dimensionless measurement strength'
    S_eff_e = None
    S_eff_g = None
    S_eff_e_selected = None
    S_eff_g_selected = None
    ################-----------------##
    ######################################################


    ######################################################
    ####   Results   ---------------##
    ###  ##########
    ### Count states ###
    dict_count = {
    'n_left_g'  : 0,    'n_left_e'  : 0,
    'n_right_g' : 0,    'n_right_e' : 0,
    'p_left_g'  : 0,    'p_left_e'  : 0,
    'p_right_g' : 0,    'p_right_e' : 0,
    'sign_ge'   : 0,
    'p_gg'      : 0,    'p_ee'      : 0,
    'p_ge'      : 0,    'p_eg'      : 0,
    }
    ### Count states after postselection ###
    dict_count_select = {
    'n_left_g'  : 0,    'n_left_e'  : 0,
    'n_right_g' : 0,    'n_right_e' : 0,
    'p_left_g'  : 0,    'p_left_e'  : 0,
    'p_right_g' : 0,    'p_right_e' : 0,
    'sign_ge'   : 0,
    'p_gg'      : 0,    'p_ee'      : 0,
    'p_ge'      : 0,    'p_eg'      : 0,
    }
    ### Fidelity dictionary ###
    dict_fidelity = None
    ################-----------------##
    ######################################################


    ############################################################################
    ############################################################################
    #### METHODS ###############################################################

    def __init__(self, data=None, datafile=None, paramfile=None, param=None, nbins=400, all_included=True):

        def get_timestamp(filename):
            '''
            function just to take string of time of measurement from file
            '''
            timestr = ''
            try:
                #opening file
                f = open(filename, "r")
                lines = f.readlines()
                f.close()
            except:
                print 'Error! get_timestamp() can not open the file'
                return None

            timestr = lines[1]
            timestr = timestr[13:-1]

            return timestr

        def get_parameters(filename):
            '''
            Takes file of parameters and extract parameters values. Returns dictionary.
            From file
            '''
            ### opening file
            try:
                f = open(filename, "r")
                lines = f.readlines()
                f.close()
            except:
                print 'Can not read the file'
                return None

            if ( len(lines)== 45 ):      ### new format of parameters file
                print '2nd generation of parameters object'
                par_dict = {'freq_read':0, 'power1':0, 't_read':0, 'rudat':0, 'freq_q':0, 'power2':0, 'rudat2':0, 'tpi':0, 'nsigma':0, 'cur':0}
                NUM_OF_VALUES = 10

                # extract numbers from string one by one
                string = lines[44]
                list_of_values_str = []
                for i in range(NUM_OF_VALUES):
                    ## it is 9 parameters to read
                    cut_end = string.find('\t')
                    str = string[ 0: cut_end ]
                    list_of_values_str.append(float(str))
                    string = string[ cut_end+1 : ]

                #fulfill dictionary and return
                par_dict['freq_read'] = list_of_values_str[0]
                par_dict['power1'] =    list_of_values_str[1]
                par_dict['t_read'] =    list_of_values_str[2]
                par_dict['rudat']  =    list_of_values_str[3]
                par_dict['freq_q'] =    list_of_values_str[4]
                par_dict['power2'] =    list_of_values_str[5]
                par_dict['rudat2'] =    list_of_values_str[6]
                par_dict['tpi']    =    list_of_values_str[7]
                par_dict['nsigma'] =    list_of_values_str[8]
                par_dict['cur']    =    list_of_values_str[9]

                return par_dict


            elif ( len(lines)== 41 ):
                print 'Caution! 1st generation of parameters object'
                par_dict = {'power1':0, 'freq_read':0, 't_read':0, 'rudat':0, 'power2':0, 'freq_q':0, 'tpi':0, 'cur':0, 'nsigma':0, 'rudat2':None}
                NUM_OF_VALUES = 9

                # extract numbers from string one by one
                string = lines[40]
                list_of_values_str = []
                for i in range(NUM_OF_VALUES):
                    ## it is 9 parameters to read
                    cut_end = string.find('\t')
                    str = string[ 0: cut_end ]
                    list_of_values_str.append(float(str))
                    string = string[ cut_end+1 : ]

                #fulfill dictionary and return
                par_dict['power1'] = list_of_values_str[0]
                par_dict['freq_read'] =    list_of_values_str[1]
                par_dict['t_read'] =    list_of_values_str[2]
                par_dict['rudat']  =    list_of_values_str[3]
                par_dict['power2'] =    list_of_values_str[4]
                par_dict['freq_q'] =    list_of_values_str[5]
                par_dict['tpi'] =    list_of_values_str[6]
                par_dict['cur']    =    list_of_values_str[7]
                par_dict['nsigma'] =    list_of_values_str[8]

                return par_dict

            elif ( len(lines)== 49):
                print '3rd generation of parameters object'
                par_dict = {'power1':0, 'freq_read':0, 't_read':0, 'rudat':0, 'power2':0, 'freq_q':0, 'tpi':0, 'cur':0, 'nsigma':0, 'rudat2':None, 'phase1':None}
                NUM_OF_VALUES = 11

                # extract numbers from string one by one
                string = lines[48]
                list_of_values_str = []
                for i in range(NUM_OF_VALUES):
                    ## it is 9 parameters to read
                    cut_end = string.find('\t')
                    str = string[ 0: cut_end ]
                    list_of_values_str.append(float(str))
                    string = string[ cut_end+1 : ]

                #fulfill dictionary and return
                par_dict['freq_read'] = list_of_values_str[0]
                par_dict['power1'] =    list_of_values_str[1]
                par_dict['t_read'] =    list_of_values_str[2]
                par_dict['rudat']  =    list_of_values_str[3]
                par_dict['freq_q'] =    list_of_values_str[4]
                par_dict['power2'] =    list_of_values_str[5]
                par_dict['rudat2'] =    list_of_values_str[6]
                par_dict['tpi']    =    list_of_values_str[7]
                par_dict['nsigma'] =    list_of_values_str[8]
                par_dict['cur']    =    list_of_values_str[9]
                par_dict['phase1']    =    list_of_values_str[10]

                return par_dict



            else:
                print 'Error of loading. Can not recognise the format of parameter file.'
                return None

        def get_paramobject(parameters):
            '''
            Takes file of parameters and extract parameters values. Returns dictionary.
            From object Parameters
            '''
            ### opening file

            # try:
            par_dict = {'freq_read':0, 'power1':0, 't_read':0, 'rudat':0, 'freq_q':0, 'power2':0, 'rudat2':0, 'tpi':0, 'nsigma':0, 'cur':0, 'phase1':0}
            NUM_OF_VALUES = 11

            #fulfill dictionary and return
            try:
                par_dict['freq_read']   = parameters.freq_read
            except:
                print 'cant load freq_read from parameters'
            try:
                par_dict['power1']      = parameters.power1
            except:
                print 'cant load power1 from parameters'
            try:
                par_dict['t_read']      = parameters.t_read
            except:
                print 'cant load t_read from parameters'
            try:
                par_dict['rudat']       = parameters.rudat
            except:
                print 'cant load rudat from parameters'
            try:
                par_dict['freq_q']      = parameters.freq_q
            except:
                print 'cant load freq_q from parameters'
            try:
                par_dict['power2']      = parameters.power2
            except:
                print 'cant load power2 from parameters'
            try:
                par_dict['rudat2']      = parameters.rudat2
            except:
                print 'cant load rudat2 from parameters'
            try:
                par_dict['tpi']         = parameters.tpi
            except:
                print 'cant load tpi from parameters'
            try:
                par_dict['nsigma']      = parameters.nsigma
            except:
                print 'cant load nsigma from parameters'
            try:
                par_dict['cur']         = parameters.current
            except:
                print 'cant load current from parameters'
            try:
                par_dict['phase1']         = parameters.phase1
            except:
                print 'cant load phase1 from parameters'

            print 'parameters loaded'

            return par_dict


        ########################################################################
        ###### ZERO ######
        ### Raw Data ###
        self.void_re = 0
        self.void_im = 0
        self.re_g     = np.array([])
        self.im_g     = np.array([])
        self.re_e     = np.array([])
        self.im_e     = np.array([])
        self.re_g_pre = np.array([])
        self.im_g_pre = np.array([])
        self.re_e_pre = np.array([])
        self.im_e_pre = np.array([])

        ### save angle of rotation and shifting the data
        ### to go from raw_basis to normal: use change_basis_point(x,y, SHIFT_X, SHIFT_Y, theta)
        ### to go from normal to raw_basis: use change_basis_point(x,y, -SHIFT_X, -SHIFT_Y, -theta)
        self.THETA = 0
        self.SHIFT_X = 0
        self.SHIFT_Y = 0

        ### Normalised Data ###
        self.void_x = 0
        self.void_y = 0
        self.x_g         = None
        self.x_e         = None
        self.x_g_pre     = None
        self.x_e_pre     = None
        self.y_g         = None
        self.y_e         = None
        self.y_g_p       = None
        self.y_e_pre     = None

        ### Normalised Data after Postselection ###
        self.x_g_select  = None
        self.x_e_select  = None
        self.y_g_select  = None
        self.y_e_select  = None

        ######################################################
        ##### CENTERS, THRESHOLD -----------------##
        ### centers of blobs (in first approximation - np.mean, after taken from gauss fit)
        center_x_g = None
        center_x_e = None
        center_x_g_select = None # it is different. I dont know why yet
        center_x_e_select = None
        sizeblob_x_g = None
        sizeblob_x_e = None
        sizeblob_y_g = None
        sizeblob_y_e = None
        squizingblob_g = None
        squizingblob_e = None
        ### crosses of size of each blolbs
        ## cross_g = [ax,ay,bx,by,cx,cy,dx,dy]
        cross_g_xy = None
        cross_e_xy = None

        # ### for future - need a back rotation #(normalisation)^-1
        center_re_g = None
        center_re_e = None
        center_im_g = None
        center_im_e = None
        ### crosses of size of each blolbs in real raw axes
        ## cross_g = [ax,ay,bx,by,cx,cy,dx,dy]
        cross_g_reim = None
        cross_e_reim = None

        ################-----------------##
        ######################################################

        ### dimensionless R-axis (like X in [V])#####-------##
        self.r_g         = None
        self.r_e         = None
        self.r_g_pre     = None
        self.r_e_pre     = None
        self.r_g_select  = None
        self.r_e_select  = None

        self.center_r_g  = None
        self.center_r_e  = None
        self.void_r      = None
        self.void_r      = None
        self.threshold_r = None

        ## Histograms ### (it is all objects of class Histograms)
        self.hist_x_g         = None
        self.hist_x_e         = None
        self.hist_x_g_pre     = None
        self.hist_x_e_pre     = None
        self.hist_x_g_select  = None
        self.hist_x_e_select  = None

        self.hist_y_g         = None
        self.hist_y_e         = None

        self.dict_fidelity = {
        'F'         : 0,
        'F_g'       : 0,
        'F_e'       : 0,
        'F_post'    : 0,
        'F_post_g'  : 0,
        'F_post_e'  : 0,
        'F_gaus'    : 0,
        'F_gaus_eg' : 0,
        'F_gaus_ge' : 0,
        'Err_e'     : 0,
        'Err_g'     : 0
        }

        self.dict_count_select = {
        'n_left_g'  : 0,    'n_left_e'  : 0,
        'n_right_g' : 0,    'n_right_e' : 0,
        'p_left_g'  : 0,    'p_left_e'  : 0,
        'p_right_g' : 0,    'p_right_e' : 0,
        'sign_ge'   : 0,
        'p_gg'      : 0,    'p_ee'      : 0,
        'p_ge'      : 0,    'p_eg'      : 0,
        }

        self.dict_count_select = {
        'n_left_g'  : 0,    'n_left_e'  : 0,
        'n_right_g' : 0,    'n_right_e' : 0,
        'p_left_g'  : 0,    'p_left_e'  : 0,
        'p_right_g' : 0,    'p_right_e' : 0,
        'sign_ge'   : 0,
        'p_gg'      : 0,    'p_ee'      : 0,
        'p_ge'      : 0,    'p_eg'      : 0,
        }

        ##############################################

        ### set up metadata
        self.datafile = datafile
        self.paramfile = paramfile
        if datafile is not None:
            self.timestamp = get_timestamp(datafile)

        ##############################################
        ##############################################

        ### set up data measurement
        if datafile is not None:
            self.load_data(datafile)
        elif data is not None:
            self.take_data(data)
        else:
            print 'To load the data you should either give a datafile=... either data= list of numpy.arrays'
            print 'Error of loading data SSResult'

        ### set up parameters
        if paramfile is not None:
            self.dict_param = get_parameters(paramfile)
        elif param is not None:
            self.dict_param = get_paramobject(param)
        else:
            self.dict_param = None
            print 'warning: no parameters is attached to SSResult object'

        ##############################################
        ### ALL THE PROCESSING IS HERE ####
        ##############################################
        if all_included:
            ### normalize it
            try:
                self.make_norm_data_from_raw()
            except:
                print '     ___SSResult error during normalisation'

            ### find the best threshold and shift the data
            try:
                self.set_best_threshold()
                self.shift_x_to_threshold_be_zero()
            except:
                print '     ___SSResult error during setting threshold'

            ### do postselection
            try:
                if self.make_postselected_data_from_norm() != False:
                ### make histograms
                    self.make_histograms(nbins = nbins)
                    self.make_histograms_y(nbins = nbins)

                ### make crosses for plot size of real blob on scattering diagram
                    self.make_blob_crosses()

                ### calculate fidelity
                    self.calculate_fidelity_post()
                else:
                    print 'smth went wrong'
            except:
                print '     ___SSResult error during postselection or hists'

            ## go from x mV to unitless variable r for S extraction
            try:
                self.make_x_dimensionless()
                self.make_unitless_histograms(nbins = nbins)
            except:
                print '     ___SSResult error during making unitless histograms'


        print 'Object is created'
    ############################################################################

    ### Process data methods ###
    def take_data(self, data):
        '''
        Load data from the measurement
        Takes list of numpy.arrays this shape
        [re_g, im_g, re_e, im_e, re_g_post, im_g_post, re_e_post, im_e_post]
        '''
        try:
            if type(data) is not list:
                print 'SSR Error: data must be a list'
                return False
            if len(data) < 4:
                print 'SSR Error: list of data is too short?'
                return False
            if len(data) > 12:
                print 'SSR Error: list of data is too long?'
                return False
            for datum in data:
                if type(datum) is not np.ndarray:
                    print 'data should be a np.array'
                    return False
        except:
            print 'error during checking!'
            print 'Error SSResult take_data()'
            return False

        try:
            ### switch from [Volts] to [mV] if it is aksed
            if self.CONVERT_TOMV == True:
                coef = 1000.0
            else:
                coef = 1.0

            # [re_g, im_g, re_e, im_e, re_g_post, im_g_post, re_e_post, im_e_post]
            self.re_g     = coef * data[0]
            self.im_g     = coef * data[1]
            self.re_e     = coef * data[2]
            self.im_e     = coef * data[3]
            self.re_g_pre = coef * data[4]
            self.im_g_pre = coef * data[5]
            self.re_e_pre = coef * data[6]
            self.im_e_pre = coef * data[7]
            return True

        except:
            print 'some error of SSResult.take_data()'
            return False

    def load_data(self, datafile):
        '''
        function open file and create data sequences
        '''
        def datatype(file):
            '''
            recognize type of datafile. Work for SingleShot files only
            '''
            try:
                f = open(file, "r")
                lines = f.readlines()
                f.close()
            except:
                print 'Error datatype() - can not open file'
                return 'error'

            if (lines[4]  == '#\tname: Real \n') and (lines[9] == '#\tname: Real-pi\n') and (lines[14] == '#\tname: Imag \n') and (lines[18] == '#\tname: Imag-pi \n'):
                return 'type1' ### old type
            elif (lines[4] == '#\tname: re_g\n') and (lines[7] == '#\tname: im_g\n') and (lines[10] == '#\tname: re_e\n') and (lines[13] == '#\tname: im_e\n'):
                return 'type2' ### new type
            else:
                return 'error'

        try:
            datatype = datatype(datafile)

            if self.CONVERT_TOMV == True:
                coef = 1000.0
            else:
                coef = 1.0

            if datatype == 'type2': ##new data
                raw_data_ss = np.loadtxt(datafile)
                self.re_g     = coef * raw_data_ss[:,0]
                self.im_g     = coef * raw_data_ss[:,1]
                self.re_e     = coef * raw_data_ss[:,2]
                self.im_e     = coef * raw_data_ss[:,3]
                self.re_g_pre = coef * raw_data_ss[:,4]
                self.im_g_pre = coef * raw_data_ss[:,5]
                self.re_e_pre = coef * raw_data_ss[:,6]
                self.im_e_pre = coef * raw_data_ss[:,7]
                print 'data loaded'
                return True
            elif  datatype == 'type1':   ##old data
                raw_data_ss = np.loadtxt(datafile)
                self.re_e     = coef * raw_data_ss[:,0]
                self.re_g     = coef * raw_data_ss[:,1]
                self.im_e     = coef * raw_data_ss[:,2]
                self.im_g     = coef * raw_data_ss[:,3]
                self.re_e_pre = coef * raw_data_ss[:,4]
                self.re_g_pre = coef * raw_data_ss[:,5]
                self.im_e_pre = coef * raw_data_ss[:,6]
                self.im_g_pre = coef * raw_data_ss[:,7]

                print 'data loaded'
                return True
            else:
                print 'ERROR! Can not recognize data type'
                return False

        except:
            print 'Warning load_data() error!'
            print 'can not load file: ', datafile
            return False

    def make_norm_data_from_raw(self):
        '''
        Make a normalisation of data on the axis between g_mean and e_mean values
        Write new data as object parameters
        Also make histograms of normalised data
        Returns True if finished
        '''
        ### check if it is data. if not - try to load from file
        if (self.re_g is not None) and (self.im_g is not None) and (self.re_e is not None) and (self.im_e is not None):
            re_g = self.re_g
            im_g = self.im_g
            re_e = self.re_e
            im_e = self.im_e
        else:
            print 'It is no raw data. \n loading...'
            success_load = self.load_data(self.datafile)
            if not success_load:
                print 'load was not successful. Sheck datafile'
                return False
            re_g = self.re_g
            im_g = self.im_g
            re_e = self.re_e
            im_e = self.im_e

        ### ----
        if self.re_g_pre is not None:
            re_g_pre = self.re_g_pre
        if self.im_g_pre is not None:
            im_g_pre = self.im_g_pre
        if self.re_e_pre is not None:
            re_e_pre = self.re_e_pre
        if self.im_e_pre is not None:
            im_e_pre = self.im_e_pre

        ##########______NORMALIZATION________###################################
        ### find centers of blobs
        [c_re_g, c_im_g, c_re_e, c_im_e ] = centers_two_blobs(re_g, im_g, re_e, im_e)

        ### find angle 2*alpha (angle between two blolbs according to void-state)
        angle_between_blobs = angle_three_points(c_re_g,c_im_g, self.void_re,self.void_im, c_re_e,c_im_e)

        ### find distance and theta between this centers
        [dist, theta] = complex_num_relationships(c_re_g,c_im_g,c_re_e,c_im_e)      #extract theta
        threshold_re = np.mean([c_re_g, c_re_e])  #x0
        threshold_im = np.mean([c_im_g, c_im_e])  #y0

        ### change the basis according to positions of blobs centers
        [ [re_g, im_g], [re_e, im_e] ]                      = change_basis_blobs_inf(threshold_re, threshold_im, theta, [re_g, im_g] , [re_e, im_e] )
        [ [re_g_pre, im_g_pre],  [re_e_pre, im_e_pre] ]     = change_basis_blobs_inf(threshold_re, threshold_im, theta, [re_g_pre, im_g_pre] , [re_e_pre, im_e_pre] )

        # normalize VOID state
        [void_re,void_im]                             = change_basis_point(self.void_re, self.void_im, threshold_re, threshold_im, theta)
        ### Calculate shift after rotation
        [x00, y00] = change_basis_point(0,0, threshold_re,threshold_im, theta)

        ########_____SAVING____________#########################################
        self.THETA   = self.THETA + theta
        self.SHIFT_X = self.SHIFT_X -x00
        self.SHIFT_Y = self.SHIFT_Y -y00
        self.void_x      = void_re
        self.void_y      = void_im
        self.x_g            = re_g
        self.x_e            = re_e
        self.x_g_pre        = re_g_pre
        self.x_e_pre        = re_e_pre
        self.y_g            = im_g
        self.y_e            = im_e
        self.y_g_pre        = im_g_pre
        self.y_e_pre        = im_e_pre

        ### new centers:
        [c_x_g, c_y_g, c_x_e, c_y_e ] = centers_two_blobs(self.x_g, self.y_g, self.x_e, self.y_e)

        self.center_x_g = c_x_g
        self.center_x_e = c_x_e
        print 'new center of blobs: ', self.center_x_g, ' ', self.center_x_e

        print 'data was normalised and saved'
        return True

    def set_best_threshold(self):
        '''
        Find best threshold,
        set it t oobject
        also calculate set self.dict_count values
        returns True if all good
        '''
        ###---
        def get_best_threshold(x_g, x_e, th0, delta0, permiss=1e-4, fragment=40, previous_fid=None):
            '''
            Searching for best value of threshold
            Bruteforce method
            permiss - permissiable value
            '''
            def get_fro_vs_threshold(x_g, x_e, threshold):
                '''
                Simplest function to calculate only f_ro for given threshold value
                Calculate fidelity for given value of THRESHOLD
                Takes 1D arrays of normalized g and e results. ( re_g & re_e )
                return only value of f_ro = 1. - 0.5*(p_ge + p_eg)
                '''
                dict_count = get_count_states(x_g,x_e,threshold)
                if dict_count is not None:
                    p_ge = dict_count['p_ge']
                    p_eg = dict_count['p_eg']
                    f_ro = 1. - 0.5*(p_ge + p_eg)
                    return f_ro
                else:
                    return 0

            # th_min = th0-delta0
            # th_max = th0+delta0
            # step = (th_max - th_min)/fragment
            # thresholds = np.arange(th_min, th_max, step)
            #
            thresholds = np.linspace(th0-delta0, th0+delta0, fragment)


            fids_list = []
            for th in thresholds:
                fid = get_fro_vs_threshold(x_g, x_e, th)
                fids_list.append(fid)

            arg_best = np.argmax(fids_list)
            best_th  = thresholds[ arg_best ]

            if (arg_best > 0):
                left_th  = thresholds[ arg_best -1 ]
            else:
                left_th  = best_th

            if (arg_best < len(thresholds)-1):
                right_th = thresholds[ arg_best +1 ]
            else:
                right_th = best_th

            best_fid  = get_fro_vs_threshold(x_g, x_e, best_th)
            left_fid  = get_fro_vs_threshold(x_g, x_e, left_th)
            right_fid = get_fro_vs_threshold(x_g, x_e, right_th)

            if previous_fid is not None:
                if abs(previous_fid - best_fid) < permiss:
                    return best_th

            # print 'arg best:', arg_best
            # print 'best fidelity: ', best_fid

            if abs(best_fid - left_fid) > permiss  or  abs(best_fid - right_fid) > permiss:
                # recursion here
                # print 'ONE MORE ITTERATION'
                best_th = get_best_threshold(x_g, x_e, best_th, delta0/10, permiss=permiss, previous_fid=best_fid)
            return best_th
        ###---
        ##########______NORMALIZE IF NECESSARY_____#############################
        if (self.x_g is None) or (self.x_e is None):
            success_norm = self.make_norm_data_from_raw()
            if not success_norm:
                print 'can not normalise. error of setting threshold'
                return False

        ##########______THRESHOLD_____##########################################
        ### Find best threshold value ###
        [leftlim, rightlim, rabbish1, rabbish2 ] = crop_fluctuations(self.x_g, [0] , self.x_e, [0], 0, 0 )
        delta = abs(rightlim - leftlim)

        threshold = get_best_threshold(self.x_g, self.x_e, 0, delta, permiss=1e-4)
        # print 'threshold:', threshold

        ##########______SET AND RETURN_____#####################################
        self.threshold = threshold
        self.dict_count = get_count_states(self.x_g,self.x_e, self.threshold)
        self.sign_ge = self.dict_count['sign_ge']

        print 'threshold set'
        return True

    def shift_x_to_threshold_be_zero(self):
        '''
        Search the best threshold value and shift the data on it
        returns threshold value, that was founded and used for shift
        '''
        if (self.x_g is None) or (self.x_e is None):
            success_norm = self.make_norm_data_from_raw()
            success_th = self.set_best_threshold()      #if we load new norm data - we need to redefine a threshold
            if (not success_norm) or (not success_th):
                print 'can not normalise or set threshold. error of shift_x_to_threshold_be_zero'
                return False

        if self.threshold is None:
            success_th = self.set_best_threshold()
            if not success_th:
                print 'can not set threshold. error of shift_x_to_threshold_be_zero'
                return False

        if self.threshold == 0:
            print 'Threshold is already zero'
            return True

        ###  NOW SHIFT THE NORMALIZED DATA FOR THRESHOLD TO BE ZERO ###

        threshold = self.threshold

        def shifter(arr, value):
            '''
            This function just shift a given array on given value
            '''
            new_arr = np.zeros_like(arr)
            for i in range(len(arr)):
                new_arr[i] = arr[i] - value
            return new_arr

        self.x_g = shifter(self.x_g, threshold)
        self.x_e = shifter(self.x_e, threshold)

        if self.x_g_pre is not None:
            self.x_g_pre = shifter(self.x_g_pre, threshold)
        if self.x_e_pre is not None:
            self.x_e_pre = shifter(self.x_e_pre, threshold)

        if self.x_g_select is not None:
            self.x_g_select = shifter(self.x_g_select, threshold)
        if self.x_e_select is not None:
            self.x_e_select = shifter(self.x_e_select, threshold)

        [self.void_x] = shifter([self.void_x], threshold)

        ### ( !V why is it commented 191023 )
        ### ( !V uncommented 191110 works same)
        # ### shift the center of blobs also
        if self.center_x_g is not None:
            [self.center_x_g] = shifter([self.center_x_g], threshold)
        if self.center_x_e is not None:
            [self.center_x_e] = shifter([self.center_x_e], threshold)
        print 'new centers of blobs:', self.center_x_g, self.center_x_e

        ### Kostil: And if the data is shifted we need also to shift histograms if it is exist!
        if (self.hist_x_g is not None) or (self.hist_x_e is not None) or (self.hist_x_g_select is not None)or (self.hist_x_e_select is not None):
            print 'we redo histograms, because of the shift'
            self.make_histograms()

        [self.SHIFT_X] = shifter([self.SHIFT_X], -threshold)

        self.threshold = 0

        print 'x-data shifted on', threshold, '. threshold=0'
        return True

    def make_x_dimensionless(self):
        '''
        Function to switch x from mv to dimensionless parameter
        It is need for calculate measurement strength and quantum efficiency
        https://arxiv.org/abs/1506.08165
        '''
        if self.center_x_g is not None and self.center_x_e is not None:
            c_g = self.center_x_g
            c_e = self.center_x_e
        else:
            print 'cant make dimensionless without centers'

        coef = 2/abs( c_g - c_e )
        shift = np.mean([c_g,c_e])

        ###----###
        ### make arrays elements unitless
        self.r_g        = (self.x_g - shift ) * coef
        self.r_e        = (self.x_e - shift ) * coef
        self.r_g_pre    = (self.x_g_pre - shift ) * coef
        self.r_e_pre    = (self.x_e_pre - shift ) * coef
        self.r_g_select = (self.x_g_select - shift ) * coef
        self.r_e_select = (self.x_e_select - shift ) * coef

        self.center_r_g = (self.center_x_g - shift ) * coef
        self.center_r_e = (self.center_x_e - shift ) * coef
        self.center_r_g_select = (self.center_x_g_select - shift ) * coef
        self.center_r_e_select = (self.center_x_e_select - shift ) * coef
        self.void_r     = (self.void_x - shift ) * coef
        self.threshold_r= (self.threshold - shift ) * coef

        return

    def make_postselected_data_from_norm(self):
        '''
        Do the postselection. Threshold usually=0 because we did set_best_threshold_as_zero() before
        also calculate and set dict_count_select
        '''
        if self.threshold is None:
            success_th = self.set_best_threshold()
            if not success_th:
                print 'threshold have not been defined yet'
                return False
        threshold = self.threshold

        if (self.x_g is None) or (self.x_e is None) or (self.x_g_pre is None) or (self.x_e_pre is None):
            success_norm = self.make_norm_data_from_raw()
            success_th = self.set_best_threshold()      #if we load new norm data - we need to redefine a threshold
            if (not success_norm) or (not success_th):
                print 'can not Postselect. error of normalisation or setting threshold'
                return False

        if ( len(self.x_g) != len(self.x_g_pre) ) or ( len(self.x_e) != len(self.x_e_pre) ):
            print 'ERROR make_postselected_data_from_norm(): Size of x_g,x_e,x_g_pre,x_e_pre is not the same  '
            return False

        ### ~~~~~~~~~~~~~~~~~~~~~~~~ ###

        x_g_selected = np.array([])
        x_e_selected = np.array([])
        index_g_wrong = []
        index_e_wrong = []

        # ## can be exchanged:
        # ind1 = np.where(Re_post > threshold)
        # Re_postselected = np.delete(Re, ind1)
        def g_state(val, threshold, sign_ge):
            if sign_ge > 0:
                if val > threshold:
                    return False
                else:
                    return True
            else:
                if val > threshold:
                    return True
                else:
                    return False

        for i in range(len(self.x_g)):
        #    if sameside(  self.x_g[i], self.x_g_pre[i], ref=threshold  ):
            if g_state(self.x_g_pre[i], self.threshold, self.sign_ge):
                x_g_selected = np.append(x_g_selected, self.x_g[i])
            else:
                index_g_wrong.append(i)

        for i in range(len(self.x_e)):
            #if not sameside(  self.x_e[i], self.x_e_pre[i], ref=threshold  ):
            if g_state(self.x_e_pre[i], self.threshold, self.sign_ge):
                x_e_selected = np.append(x_e_selected, self.x_e[i])
            else:
                index_e_wrong.append(i)

        self.x_g_select = x_g_selected
        self.x_e_select = x_e_selected
        self.dict_count_select = get_count_states(self.x_g_select, self.x_e_select, self.threshold)
        if self.dict_count_select is None:
            return False

        # return [index_g_wrong, index_e_wrong]  ### could be usefull for delete exact points from raw data also
        ### if you wanna do it -- do it right here, in this function
        ## RIGHT HERE ##

        print 'postselection done'
        return True

    def make_histograms(self, nbins=100):
        '''
        This function in charge of histograms
        of fit the gauss and plot it (remove this part)
        '''
        ########################################################################
        ### if no data - load it
        if self.x_g is None or self.x_e is None:
            success_norm = self.make_norm_data_from_raw()
            success_th = self.set_best_threshold()      #if we load new norm data - we need to redefine a threshold
            if (not success_norm) or (not success_th):
                print 'can not normalise or set threshold. error of make_histograms()'
                return False

        ### if no threshold - find one
        if self.threshold is None:
            success_th = self.set_best_threshold()
            if not success_th:
                print 'can not set threshold. error of make_histograms()'
                return False

        g_crop_s = -self.sign_ge
        e_crop_s =  self.sign_ge

        ### making hists and fit for x_g & x_e ###
        self.hist_x_g = Histogram(self.x_g, nbins = nbins)
        self.hist_x_g.fit(self.threshold, g_crop_s)
        self.hist_x_e = Histogram(self.x_e, nbins = nbins)
        self.hist_x_e.fit(self.threshold, e_crop_s)

        ### set a new center

        self.center_x_g = self.hist_x_g.gauss_param[0]
        self.center_x_e = self.hist_x_e.gauss_param[0]
        self.sizeblob_x_g = 2*self.hist_x_g.gauss_param[2]
        self.sizeblob_x_e = 2*self.hist_x_e.gauss_param[2]


        ### rotate this center to raw data scale
        [self.center_re_g, self.center_im_g] = change_basis_point(self.center_x_g,0, -self.SHIFT_X, -self.SHIFT_Y, -self.THETA)
        [self.center_re_e, self.center_im_e] = change_basis_point(self.center_x_e,0, -self.SHIFT_X, -self.SHIFT_Y, -self.THETA)


        ### making hists and fit for x_g_pre & x_e_pre ###
        if self.x_g_pre is not None:
            self.hist_x_g_pre = Histogram(self.x_g_pre, nbins = nbins)
        if self.x_e_pre is not None:
            self.hist_x_e_pre = Histogram(self.x_e_pre, nbins = nbins)

        ### making hists and fit for x_g_selected & x_e_selected ###
        if self.x_g_select is not None:
            self.hist_x_g_select = Histogram(self.x_g_select, nbins = nbins)
            self.hist_x_g_select.fit(self.threshold, g_crop_s)
            ### reset a new center
            self.center_x_g_select = self.hist_x_g_select.gauss_param[0]

        if self.x_e_select is not None:
            self.hist_x_e_select = Histogram(self.x_e_select, nbins = nbins)
            self.hist_x_e_select.fit(self.threshold, e_crop_s)
            ### reset a new center
            self.center_x_e_select = self.hist_x_e_select.gauss_param[0]


        ### making hists and fit for x_g_pre & x_e_pre ###
        print 'histograms made'
        return True

    def make_histograms_y(self, nbins=100):
        ### if no data - load it
        if self.y_g is None or self.y_e is None:
            success_norm = self.make_norm_data_from_raw()
            success_th = self.set_best_threshold()      #if we load new norm data - we need to redefine a threshold
            if (not success_norm) or (not success_th):
                print 'can not normalise or set threshold. error of make_histograms()'
                return False

        g_crop_s = -self.sign_ge
        e_crop_s =  self.sign_ge ##here it is meaningless
        ### making hists and fit for x_g & x_e ###
        self.hist_y_g = Histogram(self.y_g, nbins = nbins)
        self.hist_y_g.fit(0,0,do_crop=False)
        self.hist_y_e = Histogram(self.y_e, nbins = nbins)
        self.hist_y_e.fit(0,0,do_crop=False)

        ### set a new center
        self.center_y_g = self.hist_y_g.gauss_param[0]
        self.center_y_e = self.hist_y_e.gauss_param[0]
        self.sizeblob_y_g = 2*self.hist_y_g.gauss_param[2]
        self.sizeblob_y_e = 2*self.hist_y_e.gauss_param[2]


        ### make points for rectangle around blobs and rotate its coordinates
        ## example from make_hists()
        # [self.center_re_g, self.center_im_g] = change_basis_point(self.center_x_g,0, -self.SHIFT_X, -self.SHIFT_Y, -self.THETA)
        # [self.center_re_e, self.center_im_e] = change_basis_point(self.center_x_e,0, -self.SHIFT_X, -self.SHIFT_Y, -self.THETA)


        return True

    def make_blob_crosses(self):
        '''
        Makes a cross shape points around blob in xy_axes and reim_axes for after plot it
        Also calculate squizingblob_g and squizingblob_e
        For work it's necessary to have
        self.center_x_g
        self.center_y_g
        self.center_x_e
        self.center_y_e

        '''
        def cross_coordinztes_from_center_and_size(center_x, center_y, size_x, size_y):
            '''
            Function takes x_y_center of rectangle, and its sizes A and B
           d _____ c
            |     |
            |  +  | B
            |_____|
           a   A   b
            And convert it to coordinates of its corners (a,b,c,d): [x,y]
            '''
            ax = center_x - size_x/2
            ay = center_y
            bx = center_x + size_x/2
            by = center_y
            cx = center_x
            cy = center_y + size_y/2
            dx = center_x
            dy =center_y - size_y/2

            return [ax,ay,bx,by,cx,cy,dx,dy]

        def change_basis_cross(cross_list, shift_x, shift_y, theta):
            '''
            Change basis using change_basis_point()
            takes list of shape [x0,y0, x1,y1, x2,y2... xN,yN]
            and return it in new basis
            '''
            if len(cross_list) % 2 != 0:
                print 'Error change_basis_cross() - list must consist even number of elements'
                return np.zeros_like(cross_list)

            result_list_reim = []
            for i in np.arange(0, len(cross_list), 2):
                x = cross_list[i]
                y = cross_list[i+1]
                [re,im] = change_basis_point(x,y, shift_x, shift_y, theta)
                result_list_reim.append(re)
                result_list_reim.append(im)

            return result_list_reim

        ##[ax,ay,bx,by,cx,cy,dx,dy]
        self.cross_g_xy = cross_coordinztes_from_center_and_size(self.center_x_g, self.center_y_g, self.sizeblob_x_g, self.sizeblob_y_g)
        self.cross_e_xy = cross_coordinztes_from_center_and_size(self.center_x_e, self.center_y_e, self.sizeblob_x_e, self.sizeblob_y_e)

        try:
            ### set cross to raw basis Re-Im
            ##[aRe,aIm,bRe,bIm,cRe,cIm,dRe,dIm]
            self.cross_g_reim = change_basis_cross(self.cross_g_xy, -self.SHIFT_X, -self.SHIFT_Y, -self.THETA)
            self.cross_e_reim = change_basis_cross(self.cross_e_xy, -self.SHIFT_X, -self.SHIFT_Y, -self.THETA)
        except:
            print 'Error of setting self.cross_g_reim. Probably SHIFT or THETA didnt set'

        print 'crosses around blobs calculated...'

        try:
            ### calculate squizingblob
            self.squizingblob_g = self.sizeblob_y_g / self.sizeblob_x_g
            self.squizingblob_e = self.sizeblob_y_e / self.sizeblob_x_e
        except:
            print 'Error, squizingblob_g(e) was not calculated'

        return True

    def make_unitless_histograms(self, nbins=100):
        '''
        This function in charge of histograms
        of fit the gauss and plot it (remove this part)
        '''
        ########################################################################

        ### if no data - load it
        if (self.r_g is None) or (self.r_e is None) or (self.threshold_r is None):
            print 'can not plot unitless hists'
            return False

        g_crop_s = -self.sign_ge
        e_crop_s =  self.sign_ge

        ### making hists and fit for x_g & x_e ###
        self.hist_r_g = Histogram(self.r_g, nbins = nbins)
        self.hist_r_g.fit(self.threshold_r, g_crop_s)
        self.hist_r_e = Histogram(self.r_e, nbins = nbins)
        self.hist_r_e.fit(self.threshold_r, e_crop_s)

        # ### set a new center
        # self.center_x_g = self.hist_x_g.gauss_param[0]
        # self.center_x_e = self.hist_x_e.gauss_param[0]

        ### making hists and fit for x_g_pre & x_e_pre ###
        if self.r_g_pre is not None:
            self.hist_r_g_pre = Histogram(self.r_g_pre, nbins = nbins)
        if self.r_e_pre is not None:
            self.hist_r_e_pre = Histogram(self.r_e_pre, nbins = nbins)

        ### making hists and fit for x_g_selected & x_e_selected ###
        if self.r_g_select is not None:
            self.hist_r_g_select = Histogram(self.r_g_select, nbins = nbins)
            self.hist_r_g_select.fit(self.threshold_r, g_crop_s)
            # ### reset a new center
            # self.center_r_g_select = self.hist_x_g_select.gauss_param[0]

        if self.r_e_select is not None:
            self.hist_r_e_select = Histogram(self.r_e_select, nbins = nbins)
            self.hist_r_e_select.fit(self.threshold_r, e_crop_s)
            # ### reset a new center
            # self.center_x_e_select = self.hist_x_e_select.gauss_param[0]

        ### Record dimensionless measurement strength
        sig_e = self.hist_r_e.gauss_param[1]
        sig_g = self.hist_r_g.gauss_param[1]
        S_e = (2/sig_e)**2
        S_g = (2/sig_g)**2
        self.S_eff_e = S_e
        self.S_eff_g = S_g

        ### And for data after postselection
        sig_e_sel = self.hist_r_e_select.gauss_param[1]
        sig_g_sel = self.hist_r_g_select.gauss_param[1]
        S_e_sel = (2/sig_e_sel)**2
        S_g_sel = (2/sig_g_sel)**2
        self.S_eff_e_selected = S_e_sel
        self.S_eff_g_selected = S_g_sel

        ### making hists and fit for x_g_pre & x_e_pre ###
        print 'unitless histograms made'
        return True

    def calculate_fidelity_post(self):
        '''
        New version of fidelity calculator.
        Smart threshold by default (no fit, just bruteforce)
        '''
        ### load dictionary with p_ij raw
        if (self.dict_count['p_gg'] ==0) or (self.dict_count['p_ge'] ==0) or (self.dict_count['p_eg'] ==0) or (self.dict_count['p_ee'] ==0):
            print 'Error of calculate_fidelity_post(): no data in dict_count'
            return None
        ### count states from raw data ###
        p_gg = self.dict_count['p_gg']
        p_ge = self.dict_count['p_ge']
        p_eg = self.dict_count['p_eg']
        p_ee = self.dict_count['p_ee']


        ### load dictionary with p_ij selected
        if (self.dict_count_select['p_gg'] ==0) or (self.dict_count_select['p_ge'] ==0) or (self.dict_count_select['p_eg'] ==0) or (self.dict_count_select['p_ee'] ==0):
            print 'Error of calculate_fidelity_post(): no data in dict_count'
            return None
        ### count states with postselection ###
        p_gg_post = self.dict_count_select['p_gg']
        p_ge_post = self.dict_count_select['p_ge']
        p_eg_post = self.dict_count_select['p_eg']
        p_ee_post = self.dict_count_select['p_ee']


        ### calculate fidelities
        f_ro = 1. - 0.5*(p_ge + p_eg)
        f_g  = 1. - p_ge    #fid of preparation |g> state
        f_e  = 1. - p_eg

        f_ro_post = 1. - 0.5*(p_ge_post + p_eg_post)
        f_g_post  = 1. - p_ge_post
        f_e_post  = 1. - p_eg_post


        ## calculate gaussian overlap (using function)
        if self.hist_x_g_select is not None:
            ## use selected data if possible
            g_gaus_par = self.hist_x_g_select.gauss_param
            e_gaus_par = self.hist_x_e_select.gauss_param
        else:
            g_gaus_par = self.hist_x_g.gauss_param
            e_gaus_par = self.hist_x_e.gauss_param

        dict_overlap = get_overlap_error(g_gaus_par, e_gaus_par, self.threshold, self.sign_ge)
        f_over_g    = dict_overlap['f_o_g']
        f_over_e    = dict_overlap['f_o_e']
        f_over_tot  = dict_overlap['f_o']

        p_ge_overlap = dict_overlap['err_o_g']
        p_eg_overlap = dict_overlap['err_o_e']


        #####___SAVING_RESULT____#############
        self.dict_fidelity['F']         = f_ro
        self.dict_fidelity['F_g']       = f_g
        self.dict_fidelity['F_e']       = f_e
        self.dict_fidelity['F_post']    = f_ro_post
        self.dict_fidelity['F_post_g']  = f_g_post
        self.dict_fidelity['F_post_e']  = f_e_post
        self.dict_fidelity['F_gaus' ]   = f_over_tot
        self.dict_fidelity['F_gaus_eg'] = f_over_e
        self.dict_fidelity['F_gaus_ge'] = f_over_g
        ### strange useless parameters: (just a tradition after Remy)
        self.dict_fidelity['Err_e']     = p_ge_post - p_ge_overlap
        self.dict_fidelity['Err_g']     = p_eg_post - p_eg_overlap


        #### Return ###
        print 'Fidelity was calculated'
        # return self.dict_fidelity
        return True

    def erase_data(self, data_to_erase, selfrun=False):
        '''
        Function to get free memory.
        Erase the data of object after calculations was done
        parameter data_to_erase:
            'raw' - delete re_g, re_e, im_g, im_e and pre-selections re,im
            'norm'- delete x_g, x_e, y_g, y_e, and pre-selections x,y
            'pre' - delete only results of pre-pulses: re_g_pre, x_g_pre..
            'select' - delete selectde data
        '''
        if data_to_erase == 'raw':
            self.re_g = None
            self.re_e = None
            self.im_g = None
            self.im_e = None
            self.re_g_pre = None
            self.re_e_pre = None
            self.im_g_pre = None
            self.im_e_pre = None
            if not selfrun:
                print '__Raw data erased!'
        elif data_to_erase == 'norm':
            self.x_g = None
            self.x_e = None
            self.y_g = None
            self.y_e = None
            self.x_g_pre = None
            self.x_e_pre = None
            self.y_g_pre = None
            self.y_e_pre = None
            if not selfrun:
                print '__Normed data erased!'
        elif data_to_erase == 'pre':
            self.re_g_pre = None
            self.re_e_pre = None
            self.im_g_pre = None
            self.im_e_pre = None
            self.x_g_pre = None
            self.x_e_pre = None
            self.y_g_pre = None
            self.y_e_pre = None
            if not selfrun:
                print '__Pre-pulse data erased!'
        elif data_to_erase == 'select':
            x_g_select  = None
            x_e_select  = None
            y_g_select  = None
            y_e_select  = None
            if not selfrun:
                print '__Selected data erased!'
        elif data_to_erase == 'all':
            self.erase_data(data_to_erase='raw',   selfrun=True)
            self.erase_data(data_to_erase='norm',  selfrun=True)
            self.erase_data(data_to_erase='select',selfrun=True)
            print 'All data erased!'
        else:
            print '__Nothing erased. Check the parameter'

        return True

    ### Drawing methods ###
    # def plot_scatter_two_blob_old_color_breaken(self, norm=False, centers=None, save=False, figsize=[15,10], markersize=0.2, crosssize=10, lw=1, transpcy=100e-2,  fname='Two_blob', savepath='', fig_transp = True, show=True, limits=[None,None,None,None], crop=True, dark=True, title_str=None, font=None, zero_on_plot=False, pre_read=False):
    #     '''
    #     Plots diagramm of scattering for two blobs on the i-q plane
    #     returns limits of axis (it is used for histograms)
    #     if want to set the limits dont forget to make 'crop=False'
    #     figsize - in inches!
    #     '''
    #     if norm:
    #         if (self.x_g is None) or (self.x_e is None):
    #             self.make_norm_data_from_raw()    ### do renormalization
    #
    #         void_re = self.void_x
    #         void_im = self.void_y
    #         re_g     = self.x_g
    #         im_g     = self.y_g
    #         re_e     = self.x_e
    #         im_e     = self.y_e
    #         re_g_p = self.x_g_pre
    #         im_g_p = self.y_g_pre
    #         re_e_p = self.x_e_pre
    #         im_e_p = self.y_e_pre
    #     else:
    #         if (self.re_g is None) or (self.im_g is None) or (self.re_e is None) or (self.im_e is None):
    #             print 'It is no raw data. \n loading...'
    #             success_load = self.load_data(self.datafile)
    #             if not success_load:
    #                 print 'load was not successful'
    #                 return None
    #
    #         void_re = self.void_re
    #         void_im = self.void_im
    #         re_g     = self.re_g
    #         im_g     = self.im_g
    #         re_e     = self.re_e
    #         im_e     = self.im_e
    #         re_g_p = self.re_g_pre
    #         im_g_p = self.im_g_pre
    #         re_e_p = self.re_e_pre
    #         im_e_p = self.im_e_pre
    #         if pre_read:
    #             re_g_pre = self.re_g_pre
    #             im_g_pre = self.im_g_pre
    #             re_e_pre = self.re_e_pre
    #             im_e_pre = self.im_e_pre
    #
    #     ### setting centers
    #     if centers is not None:
    #         [c_re_g, c_im_g, c_re_e, c_im_e] = centers  ## given manually
    #     else:
    #         if not norm:
    #             if (self.center_re_g is None) or (self.center_re_e is None) or (self.center_im_g is None) or (self.center_im_e is None):
    #                 [ self.center_re_g, self.center_im_g, self.center_re_e, self.center_im_e] = centers_two_blobs(re_g, im_g, re_e, im_e)
    #             [c_re_g, c_im_g, c_re_e, c_im_e ] = [ self.center_re_g, self.center_im_g, self.center_re_e, self.center_im_e ]
    #         else:
    #             if (self.center_x_g is not None) and (self.center_x_e is not None):
    #                 [c_re_g, c_im_g, c_re_e, c_im_e ] = [ self.center_x_g, 0, self.center_x_e, 0  ]
    #             else:
    #                 print 'Error: data is not normalised'
    #
    #
    #         ### SET CUSTOM FONT ###
    #     #1
    #     # plt.rc('font', family = 'DejaVuSans')
    #     ###-----------------###
    #     #2
    #     if font is None:
    #         plt.rc('font', family = 'Verdana')
    #     # fontproperties = font
    #     ### ##############T ###
    #
    #     #### custom colors #3
    #     # vector_bw_blobs = 0.7*lw
    #     # vector_state_lw = 0.7*lw
    #     # vector_zero_lw = 0.5*lw
    #     markersize = markersize
    #     crosssize = crosssize
    #
    #     if dark:
    #         color_g = my_colors_dict['blob_g']
    #         color_e = my_colors_dict['blob_e']
    #           #None - by default
    #             ### centers of clouds
    #         # color_g_mean = '#795fd7'
    #         color_g_mean = my_colors_dict['g_state_mark']
    #         color_e_mean = my_colors_dict['e_state_mark']
    #         color_dist =    my_colors_dict['deus_ex_gold']
    #             ### vectors from void_point to centers
    #         color_g_vector = my_colors_dict['deus_ex_gold']
    #         color_e_vector = my_colors_dict['deus_ex_gold']
    #             ### zero points
    #         color_void = my_colors_dict['deus_ex_gold']
    #         color_zero = my_colors_dict['meduza_gold']
    #         color_zero_vector = my_colors_dict['meduza_gold']
    #             ### background of image
    #         fig_face_color = my_colors_dict['meduza_dark'] #this does not work
    #         fig_border_color = 'r'
    #         bg_color = 'k'
    #         grid_color =  my_colors_dict['meduza_gold']
    #         grid_transp = 0.5
    #         title_color = my_colors_dict['meduza_gold']
    #         legend_color = my_colors_dict['meduza_dark']
    #         legend_text_color = my_colors_dict['meduza_gold']
    #         legend_alpha = 0.7
    #         legend_frame_color = my_colors_dict['meduza_gold']
    #
    #         AXES_COLOR = my_colors_dict['meduza_gold']
    #         import matplotlib as mpl
    #         mpl.rc('axes', edgecolor=AXES_COLOR, labelcolor=AXES_COLOR, grid=True)
    #         mpl.rc('xtick', color=AXES_COLOR)
    #         mpl.rc('ytick', color=AXES_COLOR)
    #         mpl.rc('grid', color=AXES_COLOR)
    #     else:
    #         color_g = 'b'
    #         color_e = 'r'
    #         color_g_mean = 'midnightblue'
    #         color_e_mean = 'maroon'
    #         frame_color = 'white'
    #         bg_color = 'white'
    #         grid_color = 'lightgrey'
    #         color_dist = 'gold'
    #         color_zero = 'k'
    #         grid_transp=None
    #         color_g = 'b'
    #         color_e = 'r'
    #             ### centers of clouds
    #         # color_g_mean = '#795fd7'
    #         color_g_mean = 'midnightblue'
    #         color_e_mean = 'maroon'
    #         color_dist = 'gold'
    #             ### vectors from void_point to centers
    #         color_g_vector = 'gold'
    #         color_e_vector ='gold'
    #             ### zero points
    #         color_void = 'gold'
    #         color_zero = 'k'
    #         color_zero_vector = 'k'
    #             ### background of image
    #         fig_face_color = 'white' #this does not work
    #         fig_border_color = 'r'
    #         bg_color = 'white'
    #         grid_color =  my_colors_dict['meduza_dark']
    #         grid_transp = 0.5
    #         title_color = 'k'
    #         legend_color = 'white'
    #         legend_text_color = 'k'
    #         legend_alpha = 0.7
    #         legend_frame_color = my_colors_dict['meduza_gold']
    #
    #         import matplotlib as mpl
    #         AXES_COLOR = my_colors_dict['meduza_dark']
    #         mpl.rc('axes', edgecolor=AXES_COLOR, labelcolor=AXES_COLOR, grid=True)
    #         mpl.rc('xtick', color=AXES_COLOR)
    #         mpl.rc('ytick', color=AXES_COLOR)
    #         mpl.rc('grid', color=AXES_COLOR)
    #
    #     if save:
    #         if savepath == '':
    #             savepath='savings\\'
    #
    #     str_fidelity = ''
    #     str_params   = ''
    #     if self.dict_fidelity is not None:
    #         str_fidelity = 'F_post:'+my_stround(  100*self.dict_fidelity['F_post'],4 )+'% F:'+my_stround(  100*self.dict_fidelity['F'],4 )+'% F_gaus:'+my_stround(  100*self.dict_fidelity['F_gaus'],4 )+'%'
    #     if self.dict_param is not None:
    #         str_params = 'Rdt:'+ my_stround( self.dict_param['rudat'],4 )+'dB; t_read:'+my_stround( 1e9*self.dict_param['t_read'],3 )+'ns'
    #
    #
    #     ### find angle 2*alpha (angle between two blolbs according to void-state)
    #     angle_between_blobs = angle_three_points(c_re_g,c_im_g, void_re,void_im, c_re_e,c_im_e)
    #
    #     ## SETTING BORDERS ### (now it is outside)
    #     if not crop:
    #         [leftlim, rightlim, toplim, bottomlim ]= limits
    #         if leftlim is None:
    #             leftlim =  np.min([  0,   np.min([ re_e, re_g ])  ])
    #         if rightlim is None:
    #             rightlim = np.max([  0,   np.max([ re_e, re_g ])  ])
    #         if toplim is None:
    #             toplim =   np.max([  0,   np.max([ im_e, im_g ])  ])
    #         if bottomlim is None:
    #             bottomlim =np.min([  0,   np.min([ im_e, im_g ])  ])
    #     else:
    #         [leftlim, rightlim, toplim, bottomlim ] = crop_fluctuations(re_g, im_g ,re_e, im_e, void_re=self.void_re, void_im=self.void_im )
    #
    #     ### calculation of angle and distance between blolb centers
    #     [dist, theta] = complex_num_relationships(c_re_g,c_im_g,c_re_e,c_im_e)
    #     str_blob_place = 'Distance:'+ my_stround(dist,4) + '; Theta' +":" + my_stround(theta,2,withminus=True) + 'deg'
    #     print str_blob_place
    #
    #     ### vectors of each states from void signal
    #     [amp_g, ph_g] = complex_num_relationships(void_re, void_im, c_re_g,c_im_g)
    #     [amp_e, ph_e] = complex_num_relationships(void_re, void_im, c_re_e,c_im_e)
    #
    #     str_g_place = 'Amp:'+ my_stround(amp_g,5,withminus=True) + '; Phase:'+ my_stround(ph_g,4,withminus=True)+ 'deg'
    #     str_e_place = 'Amp:'+ my_stround(amp_e,5,withminus=True) + '; Phase:'+ my_stround(ph_e,4,withminus=True)+ 'deg'
    #     if (self.squizingblob_g is not None) and (self.squizingblob_e is not None):
    #         str_g_place = str_g_place + '; sqz:' + my_stround(self.squizingblob_g,4,withminus=False)
    #         str_e_place = str_e_place + '; sqz:' + my_stround(self.squizingblob_e,4,withminus=False)
    #
    #     amp_relation = my_stround(amp_e/amp_g,4,withminus=True)
    #     str_blob_place = str_blob_place + '; Ratio amps: '+ amp_relation
    #     print str_g_place
    #     print str_e_place
    #     print 'Ratio amps: '+ amp_relation
    #
    #     if dark:
    #         fname = fname + '_dark'
    #
    #     if figsize is None:
    #         fig = plt.figure(fname, facecolor=fig_face_color, edgecolor = fig_border_color)
    #     else:
    #         if (type(figsize) != list):
    #             print 'figsize must be a lsit'
    #             return False
    #         else:
    #             if len(figsize) != 2:
    #                 print 'figsize list must contain to numbers (x and y size)'
    #                 return False
    #         fig = plt.figure(fname, facecolor=fig_face_color, edgecolor = fig_border_color, figsize=(figsize[0],figsize[1]))
    #
    #
    #     ax = fig.add_subplot(1, 1, 1) # nrows, ncols, index
    #     ax.set_facecolor(bg_color)
    #
    #     if font is not None:
    #         plt.axis('equal',fontproperties = font)   #same step X and Y        #square axis automatic
    #     else:
    #         plt.axis('equal')
    #
    #     plt.xlim( left = leftlim, right=rightlim )
    #     plt.ylim( top =  toplim, bottom=bottomlim )
    #
    #     plt.grid(color=grid_color, alpha= grid_transp)
    #
    #     if title_str is None or title_str=='':
    #         if font is not None:
    #             plt.title(self.timestamp, color=title_color,fontproperties = font)
    #         else:
    #             plt.title(self.timestamp, color=title_color)
    #     else:
    #         if font is not None:
    #             plt.title(title_str, color=title_color,fontproperties = font)
    #         else:
    #             plt.title(title_str, color=title_color)
    #
    #     if self.CONVERT_TOMV:
    #         lab_units = '[mV]'
    #     else:
    #         lab_units = '[V]'
    #
    #     if font is not None:
    #         plt.xlabel('Re '+lab_units, fontproperties = font)
    #         plt.ylabel('Im '+lab_units, fontproperties = font)
    #     else:
    #         plt.xlabel('Re '+lab_units)
    #         plt.ylabel('Im '+lab_units)
    #
    #     plt.scatter(re_e, im_e, color=color_e, alpha=transpcy, s=markersize)
    #     plt.scatter(re_g, im_g, color=color_g, alpha=transpcy, s=markersize)
    #
    #     if pre_read:
    #         if not norm:
    #             plt.scatter(re_g_pre, im_g_pre, color='g', alpha=transpcy, s=markersize)
    #
    #     ### fake plot just for string in legend
    #     plt.plot([],[], label = str_params, visible=False)
    #     plt.plot([],[], label = str_fidelity, visible=False)
    #
    #
    #     ### real plots
    #     plt.plot([c_re_g, c_re_e], [c_im_g, c_im_e], label=str_blob_place, color=color_dist, linewidth = cs['vector_bw_blobs'])              #plot line between blobs centers
    #     if zero_on_plot:
    #         plt.plot([void_re,c_re_g], [void_im, c_im_g], color=color_g_vector, linewidth=cs['vector_state_lw'] )        #plot line from VOID to |g> blob
    #         plt.plot([void_re,c_re_e], [void_im, c_im_e], color=color_e_vector, linewidth=cs['vector_state_lw'] )        #plot line from VOID to |e> blob
    #         plt.plot([0, void_re], [0, void_im], color=color_zero_vector, linewidth=cs['vector_zero_lw'])      #plot line form [0,0] to VOID
    #
    #     plot_centers = False
    #     plt.plot([ c_re_g ],[ c_im_g ], 'X', markersize=crosssize, color=color_g_mean, label='g-state: '+str_g_place, visible=plot_centers)  #this two needs only for legend color
    #     plt.plot([ c_re_e ],[ c_im_e ], 'X', markersize=crosssize, color=color_e_mean, label='e-state: '+str_e_place, visible=plot_centers)
    #
    #     #### PLOT CROSSES AROUND EACH BLOB TO SHOW REAL SIZE ###################
    #     ###remake it with one pair of plt.plot !V !V
    #     plot_crosses = True
    #     if plot_crosses:
    #         if norm:
    #             ### g blob in XY axes
    #             if self.cross_g_xy is not None:
    #                 [ax,ay,bx,by,cx,cy,dx,dy] = self.cross_g_xy
    #                 plt.plot([ax,bx],[ay,by], color='b', lw=2)
    #                 plt.plot([cx,dx],[cy,dy], color='b', lw=2)
    #             ### e blolb in XY axes
    #             if self.cross_e_xy is not None:
    #                 [ax,ay,bx,by,cx,cy,dx,dy] = self.cross_e_xy
    #                 plt.plot([ax,bx],[ay,by], color='r', lw=2)
    #                 plt.plot([cx,dx],[cy,dy], color='r', lw=2)
    #
    #         else: ## if not nramalized
    #
    #             ### g blob in ReIm axes
    #             if self.cross_g_reim is not None:
    #                 [ax,ay,bx,by,cx,cy,dx,dy] = self.cross_g_reim
    #                 plt.plot([ax,bx],[ay,by], color='b', lw=2)
    #                 plt.plot([cx,dx],[cy,dy], color='b', lw=2)
    #             ### e blolb in ReIm axes
    #             if self.cross_e_reim is not None:
    #                 [ax,ay,bx,by,cx,cy,dx,dy] = self.cross_e_reim
    #                 plt.plot([ax,bx],[ay,by], color='r', lw=2)
    #                 plt.plot([cx,dx],[cy,dy], color='r', lw=2)
    #
    #     ########################################################################
    #
    #     if zero_on_plot:
    #         zero_label = 'Void zero: Re:'+ my_stround(void_re,5,withminus=True)+ '; Im:'+ my_stround(void_im,5,withminus=True) + '; 2*alpha='+ my_stround(angle_between_blobs,3,withminus=True)+ u"\u00b0"
    #         plt.plot([void_re],[void_im],'+', label=zero_label, color=color_void )      #coordinats of no signal VOID (global)
    #         plt.plot([0],[0],'+', color=color_zero )   #coordinats of 0 V
    #
    #     if font is not None:
    #         leg = plt.legend(fancybox=True, framealpha=legend_alpha, loc='lower left', facecolor=legend_color,edgecolor=legend_frame_color, prop=font)
    #     else:
    #         leg = plt.legend(fancybox=True, framealpha=legend_alpha, loc='lower left', facecolor=legend_color,edgecolor=legend_frame_color)
    #
    #     for text in leg.get_texts():        #set color to legend text
    #         # plt.setp(text, size = fontsize)
    #         plt.setp(text, color = legend_text_color)
    #
    #
    #     if font is not None:
    #         for label in ax.get_xticklabels():  #set font to each xtick
    #             label.set_fontproperties(font)
    #         for label in ax.get_yticklabels():  #set font to each xtick
    #             label.set_fontproperties(font)
    #
    #     if save:
    #         import os
    #
    #         if not os.path.exists(savepath):
    #             os.makedirs(savepath)
    #         full_fname = savepath +'\\'+ fname + '.png'
    #         plt.savefig(full_fname, transparent = fig_transp, facecolor=fig_face_color, edgecolor=fig_border_color, fontproperties = font)
    #
    #     if show:
    #         plt.show()
    #         return [fig, ax]
    #     else:
    #         plt.close()
    #         return True

    def plot_scatter_two_blob(self, norm=False, centers=None, save=False, figsize=[15,10], transpcy=100e-2,  fname='Two_blob', savepath='', fig_transp = True, show=True, limits=[None,None,None,None], crop=True, dark=True, title_str=None, font=None, zero_on_plot=False, pre_read=False):
        '''
        Plots diagramm of scattering for two blobs on the i-q plane
        returns limits of axis (it is used for histograms)
        if want to set the limits dont forget to make 'crop=False'
        figsize - in inches!
        '''
        if norm:
            if (self.x_g is None) or (self.x_e is None):
                self.make_norm_data_from_raw()    ### do renormalization

            void_re = self.void_x
            void_im = self.void_y
            re_g     = self.x_g
            im_g     = self.y_g
            re_e     = self.x_e
            im_e     = self.y_e
            re_g_p = self.x_g_pre
            im_g_p = self.y_g_pre
            re_e_p = self.x_e_pre
            im_e_p = self.y_e_pre
        else:
            if (self.re_g is None) or (self.im_g is None) or (self.re_e is None) or (self.im_e is None):
                print 'It is no raw data. \n loading...'
                success_load = self.load_data(self.datafile)
                if not success_load:
                    print 'load was not successful'
                    return None

            void_re = self.void_re
            void_im = self.void_im
            re_g     = self.re_g
            im_g     = self.im_g
            re_e     = self.re_e
            im_e     = self.im_e
            re_g_p = self.re_g_pre
            im_g_p = self.im_g_pre
            re_e_p = self.re_e_pre
            im_e_p = self.im_e_pre
            if pre_read:
                re_g_pre = self.re_g_pre
                im_g_pre = self.im_g_pre
                re_e_pre = self.re_e_pre
                im_e_pre = self.im_e_pre

        ### setting centers
        if centers is not None:
            [c_re_g, c_im_g, c_re_e, c_im_e] = centers  ## given manually
        else:
            if not norm:
                if (self.center_re_g is None) or (self.center_re_e is None) or (self.center_im_g is None) or (self.center_im_e is None):
                    [ self.center_re_g, self.center_im_g, self.center_re_e, self.center_im_e] = centers_two_blobs(re_g, im_g, re_e, im_e)
                [c_re_g, c_im_g, c_re_e, c_im_e ] = [ self.center_re_g, self.center_im_g, self.center_re_e, self.center_im_e ]
            else:
                if (self.center_x_g is not None) and (self.center_x_e is not None):
                    [c_re_g, c_im_g, c_re_e, c_im_e ] = [ self.center_x_g, 0, self.center_x_e, 0  ]
                else:
                    print 'Error: data is not normalised'


            ### SET CUSTOM FONT ###
        #1
        # plt.rc('font', family = 'DejaVuSans')
        ###-----------------###
        #2
        if font is None:
            plt.rc('font', family = 'Verdana')
        # fontproperties = font
        ### ##############T ###

        ### set color_scheme (cs)
        if dark:
            cs = dark_scheme
        else:
            cs = bright_scheme

        ### set colors of axes
        import matplotlib as mpl
        mpl.rc('axes', edgecolor=cs['AXES_COLOR'], labelcolor=cs['AXES_COLOR'], grid=True)
        mpl.rc('xtick', color=cs['AXES_COLOR'])
        mpl.rc('ytick', color=cs['AXES_COLOR'])
        mpl.rc('grid', color=cs['AXES_COLOR'])


        if save:
            if savepath == '':
                savepath='savings\\'

        str_fidelity = ''
        str_params   = ''
        if self.dict_fidelity is not None:
            str_fidelity = 'F_post:'+my_stround(  100*self.dict_fidelity['F_post'],4 )+'% F:'+my_stround(  100*self.dict_fidelity['F'],4 )+'% F_gaus:'+my_stround(  100*self.dict_fidelity['F_gaus'],4 )+'%'
        if self.dict_param is not None:
            str_params = 'Rdt:'+ my_stround( self.dict_param['rudat'],4 )+'dB; t_read:'+my_stround( 1e9*self.dict_param['t_read'],3 )+'ns'


        ### find angle 2*alpha (angle between two blolbs according to void-state)
        angle_between_blobs = angle_three_points(c_re_g,c_im_g, void_re,void_im, c_re_e,c_im_e)

        ## SETTING BORDERS ### (now it is outside)
        if not crop:
            [leftlim, rightlim, toplim, bottomlim ]= limits
            if leftlim is None:
                leftlim =  np.min([  0,   np.min([ re_e, re_g ])  ])
            if rightlim is None:
                rightlim = np.max([  0,   np.max([ re_e, re_g ])  ])
            if toplim is None:
                toplim =   np.max([  0,   np.max([ im_e, im_g ])  ])
            if bottomlim is None:
                bottomlim =np.min([  0,   np.min([ im_e, im_g ])  ])
        else:
            [leftlim, rightlim, toplim, bottomlim ] = crop_fluctuations(re_g, im_g ,re_e, im_e, void_re=self.void_re, void_im=self.void_im )

        ### calculation of angle and distance between blolb centers
        [dist, theta] = complex_num_relationships(c_re_g,c_im_g,c_re_e,c_im_e)
        str_blob_place = 'Distance:'+ my_stround(dist,4) + '; Theta' +":" + my_stround(theta,2,withminus=True) + 'deg'
        print str_blob_place

        ### vectors of each states from void signal
        [amp_g, ph_g] = complex_num_relationships(void_re, void_im, c_re_g,c_im_g)
        [amp_e, ph_e] = complex_num_relationships(void_re, void_im, c_re_e,c_im_e)

        str_g_place = 'Amp:'+ my_stround(amp_g,5,withminus=True) + '; Phase:'+ my_stround(ph_g,4,withminus=True)+ 'deg'
        str_e_place = 'Amp:'+ my_stround(amp_e,5,withminus=True) + '; Phase:'+ my_stround(ph_e,4,withminus=True)+ 'deg'
        if (self.squizingblob_g is not None) and (self.squizingblob_e is not None):
            str_g_place = str_g_place + '; sqz:' + my_stround(self.squizingblob_g,4,withminus=False)
            str_e_place = str_e_place + '; sqz:' + my_stround(self.squizingblob_e,4,withminus=False)

        amp_relation = my_stround(amp_e/amp_g,4,withminus=True)
        str_blob_place = str_blob_place + '; Ratio amps: '+ amp_relation
        print str_g_place
        print str_e_place
        print 'Ratio amps: '+ amp_relation

        if dark:
            fname = fname + '_dark'

        if figsize is None:
            fig = plt.figure(fname, facecolor=cs['fig_face_color'], edgecolor = cs['fig_border_color'])
        else:
            if (type(figsize) != list):
                print 'figsize must be a lsit'
                return False
            else:
                if len(figsize) != 2:
                    print 'figsize list must contain to numbers (x and y size)'
                    return False
            fig = plt.figure(fname, facecolor=cs['fig_face_color'], edgecolor = cs['fig_border_color'], figsize=(figsize[0],figsize[1]))


        ax = fig.add_subplot(1, 1, 1) # nrows, ncols, index
        ax.set_facecolor(cs['bg_color'])

        if font is not None:
            plt.axis('equal',fontproperties = font)   #same step X and Y        #square axis automatic
        else:
            plt.axis('equal')

        plt.xlim( left = leftlim, right=rightlim )
        plt.ylim( top =  toplim, bottom=bottomlim )

        plt.grid(color=cs['grid_color'], alpha= cs['grid_transp'])

        if title_str is None or title_str=='':
            if font is not None:
                plt.title(self.timestamp, color=cs['title_color'], fontproperties = font)
            else:
                plt.title(self.timestamp, color=cs['title_color'])
        else:
            if font is not None:
                plt.title(title_str, color=cs['title_color'], fontproperties = font)
            else:
                plt.title(title_str, color=cs['title_color'])

        if self.CONVERT_TOMV:
            lab_units = '[mV]'
        else:
            lab_units = '[V]'

        if font is not None:
            plt.xlabel('Re '+lab_units, fontproperties = font)
            plt.ylabel('Im '+lab_units, fontproperties = font)
        else:
            plt.xlabel('Re '+lab_units)
            plt.ylabel('Im '+lab_units)

        plt.scatter(re_e, im_e, color=cs['blob_e'], alpha=transpcy, s=cs['point_scatter_size'])
        plt.scatter(re_g, im_g, color=cs['blob_g'], alpha=transpcy, s=cs['point_scatter_size'])

        if pre_read:
            if not norm:
                plt.scatter(re_g_pre, im_g_pre, color=cs['blob_prepulse'], alpha=transpcy, s=cs['point_scatter_size'])

        ### fake plot just for string in legend
        plt.plot([],[], label = str_params, visible=False)
        plt.plot([],[], label = str_fidelity, visible=False)


        ### real plots
        plt.plot([c_re_g, c_re_e], [c_im_g, c_im_e], label=str_blob_place, color=cs['color_dist'], linewidth = cs['vector_bw_blobs'] )              #plot line between blobs centers
        if zero_on_plot:
            plt.plot([void_re,c_re_g], [void_im, c_im_g], color=cs['color_g_vector'], linewidth=cs['vector_state_lw'])        #plot line from VOID to |g> blob
            plt.plot([void_re,c_re_e], [void_im, c_im_e], color=cs['color_e_vector'], linewidth=cs['vector_state_lw'])        #plot line from VOID to |e> blob
            plt.plot([0, void_re], [0, void_im], color=cs['color_zero_vector'], linewidth=cs['vector_zero_lw'])      #plot line form [0,0] to VOID

        plot_centers = False
        plt.plot([ c_re_g ],[ c_im_g ], 'X', markersize=cs['onplot_mark_size'], color=cs['color_g_mean'], label='g-state: '+str_g_place, visible=plot_centers)  #this two needs only for legend color
        plt.plot([ c_re_e ],[ c_im_e ], 'X', markersize=cs['onplot_mark_size'], color=cs['color_e_mean'], label='e-state: '+str_e_place, visible=plot_centers)

        #### PLOT CROSSES AROUND EACH BLOB TO SHOW REAL SIZE ###################
        ###remake it with one pair of plt.plot !V !V
        plot_crosses = True
        if plot_crosses:
            if norm:
                ### g blob in XY axes
                if self.cross_g_xy is not None:
                    [ax,ay,bx,by,cx,cy,dx,dy] = self.cross_g_xy
                    plt.plot([ax,bx],[ay,by], color=cs['color_g_cross'], lw=cs['lw_cross'])
                    plt.plot([cx,dx],[cy,dy], color=cs['color_g_cross'], lw=cs['lw_cross'])
                ### e blolb in XY axes
                if self.cross_e_xy is not None:
                    [ax,ay,bx,by,cx,cy,dx,dy] = self.cross_e_xy
                    plt.plot([ax,bx],[ay,by], color=cs['color_e_cross'], lw=cs['lw_cross'])
                    plt.plot([cx,dx],[cy,dy], color=cs['color_e_cross'], lw=cs['lw_cross'])

            else: ## if not nramalized

                ### g blob in ReIm axes
                if self.cross_g_reim is not None:
                    [ax,ay,bx,by,cx,cy,dx,dy] = self.cross_g_reim
                    plt.plot([ax,bx],[ay,by], color=cs['color_g_cross'], lw=cs['lw_cross'])
                    plt.plot([cx,dx],[cy,dy], color=cs['color_g_cross'], lw=cs['lw_cross'])
                ### e blolb in ReIm axes
                if self.cross_e_reim is not None:
                    [ax,ay,bx,by,cx,cy,dx,dy] = self.cross_e_reim
                    plt.plot([ax,bx],[ay,by], color=cs['color_e_cross'], lw=cs['lw_cross'])
                    plt.plot([cx,dx],[cy,dy], color=cs['color_e_cross'], lw=cs['lw_cross'])

        ########################################################################

        if zero_on_plot:
            zero_label = 'Void zero: Re:'+ my_stround(void_re,5,withminus=True)+ '; Im:'+ my_stround(void_im,5,withminus=True) + '; 2*alpha='+ my_stround(angle_between_blobs,3,withminus=True)+ u"\u00b0"
            plt.plot([void_re],[void_im],'+', label=zero_label, color=cs['color_void'], markersize=cs['onplot_mark_size'] )      #coordinats of no signal VOID (global)
            plt.plot([0],[0],'+', color=cs['color_zero'], markersize=cs['onplot_mark_size'] )   #coordinats of 0 V

        if font is not None:
            leg = plt.legend(fancybox=True, framealpha=cs['legend_alpha'], loc='lower left', facecolor=cs['legend_color'], edgecolor=cs['legend_frame_color'], prop=font)
        else:
            leg = plt.legend(fancybox=True, framealpha=cs['legend_alpha'], loc='lower left', facecolor=cs['legend_color'],edgecolor=cs['legend_frame_color'])

        for text in leg.get_texts():        #set color to legend text
            # plt.setp(text, size = fontsize)
            plt.setp(text, color = cs['legend_text_color'])


        if font is not None:
            for label in ax.get_xticklabels():  #set font to each xtick
                label.set_fontproperties(font)
            for label in ax.get_yticklabels():  #set font to each xtick
                label.set_fontproperties(font)

        if save:
            import os

            if not os.path.exists(savepath):
                os.makedirs(savepath)
            full_fname = savepath +'\\'+ fname + '.png'
            plt.savefig(full_fname, transparent = fig_transp, facecolor=cs['fig_face_color'], edgecolor=cs['fig_border_color'], fontproperties = font)

        if show:
            plt.show()
            return [fig, ax]
        else:
            plt.close()
            return True


    def plot_hists(self, regime='raw_data', dark=True, log=True, save=False, figsize=[15,10], savepath='', fname='Hists', lw=1, fig_transp=False, title_str='', font=None, show=True, limits=[None,None]):
        '''
        function of plot histograms of object is it exists
        have different regimes:
        regime = 'raw_data' - plot hists of data with fit(if it is)
        regime = 'raw_data_and_pre' - same, but with results of prepulse
        regime = 'selected' - plot postselected data with fit(if it is)
        regime = 'raw_and_selected' compare between raw data and postselected (no fit)
        '''
        ###############################################
        ### all design here _________________________##
        ###############################################

        ### STYLE
        if dark:
            cs = dark_scheme
        else:
            cs = bright_scheme

        if dark:
            fname = fname + '_dark'
            color_g = my_colors_dict['blob_g']
            color_e = my_colors_dict['blob_e']
            color_post = my_colors_dict['blob_post']
            color_fit_g = my_colors_dict['g_state_mark']
            color_fit_e = my_colors_dict['e_state_mark']
            color_g_mean = my_colors_dict['g_state_mark']
            color_e_mean = my_colors_dict['e_state_mark']

            # color_dist =    my_colors_dict['deus_ex_gold']
            # color_zero = my_colors_dict['meduza_gold']
                    ### background of image
            fig_face_color = my_colors_dict['meduza_dark'] #this does not work
            fig_border_color = 'r'
            # bg_color = 'k'
            bg_color = my_colors_dict['meduza_dark']
            grid_color =  my_colors_dict['meduza_gold']
            grid_transp = 0.5
            title_color = my_colors_dict['meduza_gold']
            legend_color = my_colors_dict['meduza_dark']
            legend_text_color = my_colors_dict['meduza_gold']
            legend_alpha = 0.7
            legend_frame_color = my_colors_dict['meduza_gold']

            th_alpha = 0.7
            th_color = my_colors_dict['deus_ex_gold']
            th_width = 2.0*lw
            th_linestyle = '--'

            AXES_COLOR = my_colors_dict['meduza_gold']
            import matplotlib as mpl
            mpl.rc('axes', edgecolor=AXES_COLOR, labelcolor=AXES_COLOR, grid=True)
            mpl.rc('xtick', color=AXES_COLOR)
            mpl.rc('ytick', color=AXES_COLOR)
            mpl.rc('grid', color=AXES_COLOR)
        else:
            color_g = 'b'
            color_e = 'r'
            color_g_mean = 'midnightblue'
            color_e_mean = 'maroon'
            color_post = my_colors_dict['blob_post']
            color_fit_g = my_colors_dict['g_state_mark']
            color_fit_e = my_colors_dict['e_state_mark']

            frame_color = 'white'
            bg_color = 'white'
            grid_color = 'lightgrey'
            color_dist = 'gold'
            color_zero = 'k'
            grid_transp=None

            color_g = 'b'
            color_e = 'r'
            transpcy=2.5e-2
            markersize = None   #None - by default
                ### centers of clouds
            # color_g_mean = '#795fd7'
            color_g_mean = 'midnightblue'
            color_e_mean = 'maroon'
            color_dist = 'gold'
            vector_bw_blobs = 0.7
                ### vectors from void_point to centers
            color_g_vector = 'gold'
            color_e_vector ='gold'
            vector_state_lw = 0.7
                ### zero points
            color_void = 'gold'
            color_zero = 'k'
            color_zero_vector = 'k'
            vector_zero_lw = 0.5
                ### background of image
            fig_face_color = 'white' #this does not work
            fig_border_color = 'r'
            bg_color = 'white'
            grid_color =  my_colors_dict['meduza_dark']
            grid_transp = 0.5
            title_color = 'k'
            legend_color = 'white'
            legend_text_color = 'k'
            legend_alpha = 0.7
            legend_frame_color = 'k'

            th_alpha = 0.7
            th_color = my_colors_dict['deus_ex_gold']
            th_width = 2.0
            th_linestyle = '--'


            import matplotlib as mpl
            AXES_COLOR = my_colors_dict['meduza_dark']
            mpl.rc('axes', edgecolor=AXES_COLOR, labelcolor=AXES_COLOR, grid=True)
            mpl.rc('xtick', color=AXES_COLOR)
            mpl.rc('ytick', color=AXES_COLOR)
            mpl.rc('grid', color=AXES_COLOR)

        ### font default
        if font is None:
            plt.rc('font', family = 'Verdana')

        ###############################################
        ### all data plot here ______________________##
        ###############################################
        if figsize is None:
            fig, ax = plt.subplots(1, 1, sharey=True, tight_layout=True, facecolor=fig_face_color, edgecolor = fig_border_color)
        else:
            if (type(figsize) != list):
                print 'figsize must be a lsit'
                return False
            else:
                if len(figsize) != 2:
                    print 'figsize list must contain to numbers (x and y size)'
                    return False
            fig, ax = plt.subplots(1, 1, figsize=(figsize[0],figsize[1]), sharey=True, tight_layout=True, facecolor=fig_face_color, edgecolor = fig_border_color)



        maxval = 1 ##this variable we use for define plt.ylim()

        if self.threshold is not None:
            plt.axvline(x=self.threshold, alpha=th_alpha, c=th_color, lw=th_width, ls=th_linestyle)

        if regime == 'raw_data':
            print 'regime: raw_data'
            if (self.center_x_g is not None) and (self.center_x_e is not None):
                plt.axvline(x=self.center_x_g, alpha=th_alpha, c=color_g_mean, lw=th_width, ls=th_linestyle)
                plt.axvline(x=self.center_x_e, alpha=th_alpha, c=color_e_mean, lw=th_width, ls=th_linestyle)
            if (self.sizeblob_x_g is not None) and (self.sizeblob_x_e is not None):
                plt.axvline(x=self.center_x_g - self.sizeblob_x_g/2, alpha=th_alpha/2, c=color_g_mean, lw=th_width, ls=th_linestyle)
                plt.axvline(x=self.center_x_g + self.sizeblob_x_g/2, alpha=th_alpha/2, c=color_g_mean, lw=th_width, ls=th_linestyle)
                plt.axvline(x=self.center_x_e - self.sizeblob_x_e/2, alpha=th_alpha/2, c=color_e_mean, lw=th_width, ls=th_linestyle)
                plt.axvline(x=self.center_x_e + self.sizeblob_x_e/2, alpha=th_alpha/2, c=color_e_mean, lw=th_width, ls=th_linestyle)

            #### plot x_g, x_e hists ####
            if (self.hist_x_g.hist_xy is not None) and (self.hist_x_e.hist_xy is not None):
                hist_g   = self.hist_x_g.hist_xy
                hist_e   = self.hist_x_e.hist_xy
                maxval = np.max([ maxval, np.max(hist_g[0]), np.max(hist_e[0]) ])
                plt.plot(hist_g[1], hist_g[0], drawstyle='steps', lw=1*lw, color=color_g, label='Read g-state')
                plt.plot(hist_e[1], hist_e[0], drawstyle='steps', lw=1*lw, color=color_e, label='Read e-state')

            #### plot fit hists ####
            if (self.hist_x_g.gauss_fit is not None) and (self.hist_x_e.gauss_fit is not None):
                hist_g  = self.hist_x_g.gauss_fit
                hist_e  = self.hist_x_e.gauss_fit
                # maxval = np.max([ maxval, np.max(hist_g_fit[0]), np.max(hist_e_fit[0]) ])
                plt.plot(hist_g[1], hist_g[0], lw=1*lw, color=color_fit_g)
                plt.plot(hist_e[1], hist_e[0], lw=1*lw, color=color_fit_e)

        elif regime == 'raw_data_and_pre':
            print 'regime: raw_data_and_pre'
            if (self.center_x_g is not None) and (self.center_x_e is not None):
                plt.axvline(x=self.center_x_g, alpha=th_alpha, c=color_g_mean, lw=th_width, ls=th_linestyle)
                plt.axvline(x=self.center_x_e, alpha=th_alpha, c=color_e_mean, lw=th_width, ls=th_linestyle)
            #### plot PREPULSE hists ####
            if (self.hist_x_g_pre is not None) and (self.hist_x_e_pre is not None):
                hist_g  = self.hist_x_g_pre.hist_xy
                hist_e  = self.hist_x_e_pre.hist_xy
                # maxval = np.max([ maxval, np.max(hist_g_fit[0]), np.max(hist_e_fit[0]) ])
                plt.plot(hist_g[1], hist_g[0], drawstyle='steps', lw=1*lw, color=color_post, label='Prepulse g-state')
                plt.plot(hist_e[1], hist_e[0], drawstyle='steps', lw=1*lw, color=color_post, label='Prepulse e-state')

            #### plot x_g, x_e hists ####
            if (self.hist_x_g.hist_xy is not None) and (self.hist_x_e.hist_xy is not None):
                hist_g   = self.hist_x_g.hist_xy
                hist_e   = self.hist_x_e.hist_xy
                maxval = np.max([ maxval, np.max(hist_g[0]), np.max(hist_e[0]) ])
                plt.plot(hist_g[1], hist_g[0], drawstyle='steps', lw=1*lw, color=color_g, label='Read g-state')
                plt.plot(hist_e[1], hist_e[0], drawstyle='steps', lw=1*lw, color=color_e, label='Read e-state')

            #### plot fit hists ####
            if (self.hist_x_g.gauss_fit is not None) and (self.hist_x_e.gauss_fit is not None):
                hist_g  = self.hist_x_g.gauss_fit
                hist_e  = self.hist_x_e.gauss_fit
                # maxval = np.max([ maxval, np.max(hist_g_fit[0]), np.max(hist_e_fit[0]) ])
                plt.plot(hist_g[1], hist_g[0], lw=1*lw, color=color_fit_g)
                plt.plot(hist_e[1], hist_e[0], lw=1*lw, color=color_fit_e)


        elif regime=='selected':
            if (self.center_x_g_select is not None) and (self.center_x_e_select is not None):
                plt.axvline(x=self.center_x_g_select, alpha=th_alpha, c=color_g_mean, lw=th_width, ls=th_linestyle)
                plt.axvline(x=self.center_x_e_select, alpha=th_alpha, c=color_e_mean, lw=th_width, ls=th_linestyle)
            if (self.sizeblob_x_g is not None) and (self.sizeblob_x_e is not None):
                plt.axvline(x=self.center_x_g - self.sizeblob_x_g/2, alpha=th_alpha/2, c=color_g_mean, lw=th_width, ls=th_linestyle)
                plt.axvline(x=self.center_x_g + self.sizeblob_x_g/2, alpha=th_alpha/2, c=color_g_mean, lw=th_width, ls=th_linestyle)
                plt.axvline(x=self.center_x_e - self.sizeblob_x_e/2, alpha=th_alpha/2, c=color_e_mean, lw=th_width, ls=th_linestyle)
                plt.axvline(x=self.center_x_e + self.sizeblob_x_e/2, alpha=th_alpha/2, c=color_e_mean, lw=th_width, ls=th_linestyle)

            print 'regime: selected'
            #### plot x_g, x_e hists  SELECTED ####
            if (self.hist_x_g_select.hist_xy is not None) and (self.hist_x_e_select.hist_xy is not None):
                hist_g   = self.hist_x_g_select.hist_xy
                hist_e   = self.hist_x_e_select.hist_xy
                maxval = np.max([ maxval, np.max(hist_g[0]), np.max(hist_e[0]) ])
                plt.plot(hist_g[1], hist_g[0], drawstyle='steps', lw=1*lw, color=color_g, label='Selected g-state')
                plt.plot(hist_e[1], hist_e[0], drawstyle='steps', lw=1*lw, color=color_e, label='Selected e-state')

            #### plot fit hists ####
            if (self.hist_x_g_select.gauss_fit is not None) and (self.hist_x_e_select.gauss_fit is not None):
                hist_g  = self.hist_x_g_select.gauss_fit
                hist_e  = self.hist_x_e_select.gauss_fit
                # maxval = np.max([ maxval, np.max(hist_g_fit[0]), np.max(hist_e_fit[0]) ])
                plt.plot(hist_g[1], hist_g[0], lw=1*lw, color=color_fit_g, label='Fit g-state selected')
                plt.plot(hist_e[1], hist_e[0], lw=1*lw, color=color_fit_e, label='Fit e-state selected')


        elif regime=='raw_and_selected':
            if (self.center_x_g_select is not None) and (self.center_x_e_select is not None):
                plt.axvline(x=self.center_x_g_select, alpha=th_alpha, c=color_g_mean, lw=th_width, ls=th_linestyle)
                plt.axvline(x=self.center_x_e_select, alpha=th_alpha, c=color_e_mean, lw=th_width, ls=th_linestyle)
            print 'regime==raw_and_selected'
            if (self.sizeblob_x_g is not None) and (self.sizeblob_x_e is not None):
                plt.axvline(x=self.center_x_g - self.sizeblob_x_g/2, alpha=th_alpha/2, c=color_g_mean, lw=th_width, ls=th_linestyle)
                plt.axvline(x=self.center_x_g + self.sizeblob_x_g/2, alpha=th_alpha/2, c=color_g_mean, lw=th_width, ls=th_linestyle)
                plt.axvline(x=self.center_x_e - self.sizeblob_x_e/2, alpha=th_alpha/2, c=color_e_mean, lw=th_width, ls=th_linestyle)
                plt.axvline(x=self.center_x_e + self.sizeblob_x_e/2, alpha=th_alpha/2, c=color_e_mean, lw=th_width, ls=th_linestyle)

            #### plot x_g, x_e hists  SELECTED ####
            if (self.hist_x_g_select.hist_xy is not None) and (self.hist_x_e_select.hist_xy is not None):
                hist_g   = self.hist_x_g_select.hist_xy
                hist_e   = self.hist_x_e_select.hist_xy
                maxval = np.max([ maxval, np.max(hist_g[0]), np.max(hist_e[0]) ])
                plt.plot(hist_g[1], hist_g[0], drawstyle='steps', lw=1*lw, color='b', label='Selected g-state')
                plt.plot(hist_e[1], hist_e[0], drawstyle='steps', lw=1*lw, color='r', label='Selected e-state')

            #### plot x_g, x_e hists ####
            if (self.hist_x_g.hist_xy is not None) and (self.hist_x_e.hist_xy is not None):
                hist_g   = self.hist_x_g.hist_xy
                hist_e   = self.hist_x_e.hist_xy
                maxval = np.max([ maxval, np.max(hist_g[0]), np.max(hist_e[0]) ])
                plt.plot(hist_g[1], hist_g[0], drawstyle='steps', lw=1*lw, color=color_g, label='Read g-state')
                plt.plot(hist_e[1], hist_e[0], drawstyle='steps', lw=1*lw, color=color_e, label='Read e-state')


        else:
            print 'plot_hist: wrong value for regime!'

        ###############################################


        ###############################################
        ### all work with axes here _________________##
        ###############################################
        ax.set_facecolor(bg_color)
        plt.grid(color=grid_color, alpha= grid_transp)

        [leftlim, rightlim] = limits
        if leftlim is None:
            leftlim = np.min([ np.min(hist_g[1]), np.min(hist_e[1]) ])
        if rightlim is None:
            rightlim = np.max([ np.max(hist_g[1]), np.max(hist_e[1]) ])
        plt.xlim(leftlim,rightlim)

        if log:
            plt.yscale('log')
            plt.ylim(ymin=1, ymax=maxval*1.5)
        else:
            plt.yscale('linear')
            plt.ylim(ymin=1, ymax=maxval*1.1)

        if title_str is '':
            title_str = self.timestamp

        if self.CONVERT_TOMV:
            lab_units = '[mV]'
        else:
            lab_units = '[V]'

        if font is not None:
            plt.title(title_str, color=title_color,fontproperties = font)
            plt.xlabel(lab_units, fontproperties = font)
            plt.ylabel('Counts',fontproperties = font)
            leg = plt.legend(fancybox=True, framealpha=legend_alpha, loc='upper left', facecolor=legend_color, edgecolor=legend_frame_color,prop=font)
            for label in ax.get_xticklabels():  #set font to each xtick
                label.set_fontproperties(font)
            for label in ax.get_yticklabels():  #set font to each xtick
                label.set_fontproperties(font)
        else:
            plt.title(title_str, color=title_color)
            plt.xlabel(lab_units)
            plt.ylabel('Counts')
            leg = plt.legend(fancybox=True, framealpha=legend_alpha, loc='upper left', facecolor=legend_color, edgecolor=legend_frame_color)

        for text in leg.get_texts():        #set color to legend text
            plt.setp(text, color = legend_text_color)

        if save:
            if savepath == '':
                savepath='savings\\'
            import os
            if not os.path.exists(savepath):
                os.makedirs(savepath)
            full_fname = savepath +'\\'+ fname + '.png'
            if fig_transp:
                plt.savefig(full_fname, transparent = True)
            else:
                plt.savefig(full_fname,facecolor=fig_face_color, edgecolor=fig_border_color)

        if show:
            plt.show()
            return fig
        else:
            plt.close()
            return True

    def plot_hists_y(self, dark=True, log=True, save=False, figsize=[15,10], savepath='', fname='Hists_Y', lw=1, fig_transp=False, title_str='', font=None, show=True, limits=[None,None]):
        '''
        function of plot histograms of object is it exists
        have different regimes:
        regime = 'raw_data' - plot hists of data with fit(if it is)

        '''
        ###############################################
        ### all design here _________________________##
        ###############################################

        ### STYLE
        if dark:
            fname = fname + '_dark'
            color_g = my_colors_dict['blob_g']
            color_e = my_colors_dict['blob_e']
            color_post = my_colors_dict['blob_post']
            color_fit_g = my_colors_dict['g_state_mark']
            color_fit_e = my_colors_dict['e_state_mark']
            color_g_mean = my_colors_dict['g_state_mark']
            color_e_mean = my_colors_dict['e_state_mark']

            # color_dist =    my_colors_dict['deus_ex_gold']
            # color_zero = my_colors_dict['meduza_gold']
                    ### background of image
            fig_face_color = my_colors_dict['meduza_dark'] #this does not work
            fig_border_color = 'r'
            # bg_color = 'k'
            bg_color = my_colors_dict['meduza_dark']
            grid_color =  my_colors_dict['meduza_gold']
            grid_transp = 0.5
            title_color = my_colors_dict['meduza_gold']
            legend_color = my_colors_dict['meduza_dark']
            legend_text_color = my_colors_dict['meduza_gold']
            legend_alpha = 0.7
            legend_frame_color = my_colors_dict['meduza_gold']

            th_alpha = 0.7
            th_color = my_colors_dict['deus_ex_gold']
            th_width = 2.0*lw
            th_linestyle = '--'

            AXES_COLOR = my_colors_dict['meduza_gold']
            import matplotlib as mpl
            mpl.rc('axes', edgecolor=AXES_COLOR, labelcolor=AXES_COLOR, grid=True)
            mpl.rc('xtick', color=AXES_COLOR)
            mpl.rc('ytick', color=AXES_COLOR)
            mpl.rc('grid', color=AXES_COLOR)
        else:
            color_g = 'b'
            color_e = 'r'
            color_g_mean = 'midnightblue'
            color_e_mean = 'maroon'
            color_post = my_colors_dict['blob_post']
            color_fit_g = my_colors_dict['g_state_mark']
            color_fit_e = my_colors_dict['e_state_mark']

            frame_color = 'white'
            bg_color = 'white'
            grid_color = 'lightgrey'
            color_dist = 'gold'
            color_zero = 'k'
            grid_transp=None

            color_g = 'b'
            color_e = 'r'
            transpcy=2.5e-2
            markersize = None   #None - by default
                ### centers of clouds
            # color_g_mean = '#795fd7'
            color_g_mean = 'midnightblue'
            color_e_mean = 'maroon'
            color_dist = 'gold'
            vector_bw_blobs = 0.7
                ### vectors from void_point to centers
            color_g_vector = 'gold'
            color_e_vector ='gold'
            vector_state_lw = 0.7
                ### zero points
            color_void = 'gold'
            color_zero = 'k'
            color_zero_vector = 'k'
            vector_zero_lw = 0.5
                ### background of image
            fig_face_color = 'white' #this does not work
            fig_border_color = 'r'
            bg_color = 'white'
            grid_color =  my_colors_dict['meduza_dark']
            grid_transp = 0.5
            title_color = 'k'
            legend_color = 'white'
            legend_text_color = 'k'
            legend_alpha = 0.7
            legend_frame_color = 'k'

            th_alpha = 0.7
            th_color = my_colors_dict['deus_ex_gold']
            th_width = 2.0
            th_linestyle = '--'


            import matplotlib as mpl
            AXES_COLOR = my_colors_dict['meduza_dark']
            mpl.rc('axes', edgecolor=AXES_COLOR, labelcolor=AXES_COLOR, grid=True)
            mpl.rc('xtick', color=AXES_COLOR)
            mpl.rc('ytick', color=AXES_COLOR)
            mpl.rc('grid', color=AXES_COLOR)

        ### font default
        if font is None:
            plt.rc('font', family = 'Verdana')

        ###############################################
        ### all data plot here ______________________##
        ###############################################
        if figsize is None:
            fig, ax = plt.subplots(1, 1, sharey=True, tight_layout=True, facecolor=fig_face_color, edgecolor = fig_border_color)
        else:
            if (type(figsize) != list):
                print 'figsize must be a lsit'
                return False
            else:
                if len(figsize) != 2:
                    print 'figsize list must contain to numbers (x and y size)'
                    return False
            fig, ax = plt.subplots(1, 1, figsize=(figsize[0],figsize[1]), sharey=True, tight_layout=True, facecolor=fig_face_color, edgecolor = fig_border_color)


        maxval = 1 ##this variable we use for define plt.ylim()

        if self.threshold is not None:
            plt.axvline(x=self.threshold, alpha=th_alpha, c=th_color, lw=th_width, ls=th_linestyle)

        if (self.center_y_g is not None) and (self.center_y_e is not None):
            plt.axvline(x=self.center_y_g, alpha=th_alpha, c=color_g_mean, lw=th_width, ls=th_linestyle)
            plt.axvline(x=self.center_y_e, alpha=th_alpha, c=color_e_mean, lw=th_width, ls=th_linestyle)
        if (self.sizeblob_y_g is not None) and (self.sizeblob_y_e is not None):
            plt.axvline(x=self.center_y_g - self.sizeblob_y_g/2, alpha=th_alpha/2, c=color_g_mean, lw=th_width, ls=th_linestyle)
            plt.axvline(x=self.center_y_g + self.sizeblob_y_g/2, alpha=th_alpha/2, c=color_g_mean, lw=th_width, ls=th_linestyle)
            plt.axvline(x=self.center_y_e - self.sizeblob_y_e/2, alpha=th_alpha/2, c=color_e_mean, lw=th_width, ls=th_linestyle)
            plt.axvline(x=self.center_y_e + self.sizeblob_y_e/2, alpha=th_alpha/2, c=color_e_mean, lw=th_width, ls=th_linestyle)

        #### plot x_g, x_e hists ####
        if (self.hist_y_g.hist_xy is not None) and (self.hist_y_e.hist_xy is not None):
            hist_g   = self.hist_y_g.hist_xy
            hist_e   = self.hist_y_e.hist_xy
            maxval = np.max([ maxval, np.max(hist_g[0]), np.max(hist_e[0]) ])
            plt.plot(hist_g[1], hist_g[0], drawstyle='steps', lw=1*lw, color=color_g, label='Read g-state')
            plt.plot(hist_e[1], hist_e[0], drawstyle='steps', lw=1*lw, color=color_e, label='Read e-state')

        #### plot fit hists ####
        if (self.hist_y_g.gauss_fit is not None) and (self.hist_y_e.gauss_fit is not None):
            hist_g  = self.hist_y_g.gauss_fit
            hist_e  = self.hist_y_e.gauss_fit
            # maxval = np.max([ maxval, np.max(hist_g_fit[0]), np.max(hist_e_fit[0]) ])
            plt.plot(hist_g[1], hist_g[0], lw=1*lw, color=color_fit_g)
            plt.plot(hist_e[1], hist_e[0], lw=1*lw, color=color_fit_e)

        ###############################################


        ###############################################
        ### all work with axes here _________________##
        ###############################################
        ax.set_facecolor(bg_color)
        plt.grid(color=grid_color, alpha= grid_transp)

        [leftlim, rightlim] = limits
        if leftlim is None:
            leftlim = np.min([ np.min(hist_g[1]), np.min(hist_e[1]) ])
        if rightlim is None:
            rightlim = np.max([ np.max(hist_g[1]), np.max(hist_e[1]) ])
        plt.xlim(leftlim,rightlim)

        if log:
            plt.yscale('log')
            plt.ylim(ymin=1, ymax=maxval*1.5)
        else:
            plt.yscale('linear')
            plt.ylim(ymin=1, ymax=maxval*1.1)

        if title_str is '':
            title_str = self.timestamp

        if self.CONVERT_TOMV:
            lab_units = '[mV]'
        else:
            lab_units = '[V]'

        if font is not None:
            plt.title(title_str, color=title_color,fontproperties = font)
            plt.xlabel(lab_units, fontproperties = font)
            plt.ylabel('Counts',fontproperties = font)
            leg = plt.legend(fancybox=True, framealpha=legend_alpha, loc='upper left', facecolor=legend_color, edgecolor=legend_frame_color,prop=font)
            for label in ax.get_xticklabels():  #set font to each xtick
                label.set_fontproperties(font)
            for label in ax.get_yticklabels():  #set font to each xtick
                label.set_fontproperties(font)
        else:
            plt.title(title_str, color=title_color)
            plt.xlabel(lab_units)
            plt.ylabel('Counts')
            leg = plt.legend(fancybox=True, framealpha=legend_alpha, loc='upper left', facecolor=legend_color, edgecolor=legend_frame_color)

        for text in leg.get_texts():        #set color to legend text
            plt.setp(text, color = legend_text_color)

        if save:
            if savepath == '':
                savepath='savings\\'
            import os
            if not os.path.exists(savepath):
                os.makedirs(savepath)
            full_fname = savepath +'\\'+ fname + '.png'
            if fig_transp:
                plt.savefig(full_fname, transparent = True)
            else:
                plt.savefig(full_fname,facecolor=fig_face_color, edgecolor=fig_border_color)

        if show:
            plt.show()
            return fig
        else:
            plt.close()
            return True


    def plot_hists_unitless(self, regime='raw_and_selected', dark=True, log=True, save=False, figsize=[15,10], savepath='', fname='Hists_unitless', lw=1, fig_transp=False, title_str='', font=None, show=True, limits=[None,None]):
        '''
        function of plot histograms of object is it exists
        have different regimes:
        regime = 'raw_data' - plot hists of data with fit(if it is)
        regime = 'raw_data_and_pre' - same, but with results of prepulse
        regime = 'selected' - plot postselected data with fit(if it is)
        regime = 'raw_and_selected' compare between raw data and postselected (no fit)
        '''
        ###############################################
        ### all design here _________________________##
        ###############################################

        ### STYLE
        if dark:
            fname = fname + '_dark'
            color_g = my_colors_dict['blob_g']
            color_e = my_colors_dict['blob_e']
            color_post = my_colors_dict['blob_post']
            color_fit_g = my_colors_dict['g_state_mark']
            color_fit_e = my_colors_dict['e_state_mark']
            color_g_mean = my_colors_dict['g_state_mark']
            color_e_mean = my_colors_dict['e_state_mark']

            # color_dist =    my_colors_dict['deus_ex_gold']
            # color_zero = my_colors_dict['meduza_gold']
                    ### background of image
            fig_face_color = my_colors_dict['meduza_dark'] #this does not work
            fig_border_color = 'r'
            # bg_color = 'k'
            bg_color = my_colors_dict['meduza_dark']
            grid_color =  my_colors_dict['meduza_gold']
            grid_transp = 0.5
            title_color = my_colors_dict['meduza_gold']
            legend_color = my_colors_dict['meduza_dark']
            legend_text_color = my_colors_dict['meduza_gold']
            legend_alpha = 0.7
            legend_frame_color = my_colors_dict['meduza_gold']

            th_alpha = 0.7
            th_color = my_colors_dict['deus_ex_gold']
            th_width = 2.0*lw
            th_linestyle = '--'

            AXES_COLOR = my_colors_dict['meduza_gold']
            import matplotlib as mpl
            mpl.rc('axes', edgecolor=AXES_COLOR, labelcolor=AXES_COLOR, grid=True)
            mpl.rc('xtick', color=AXES_COLOR)
            mpl.rc('ytick', color=AXES_COLOR)
            mpl.rc('grid', color=AXES_COLOR)
        else:
            color_g = 'b'
            color_e = 'r'
            color_g_mean = 'midnightblue'
            color_e_mean = 'maroon'
            color_post = my_colors_dict['blob_post']
            color_fit_g = my_colors_dict['g_state_mark']
            color_fit_e = my_colors_dict['e_state_mark']

            frame_color = 'white'
            bg_color = 'white'
            grid_color = 'lightgrey'
            color_dist = 'gold'
            color_zero = 'k'
            grid_transp=None

            color_g = 'b'
            color_e = 'r'
            transpcy=2.5e-2
            markersize = None   #None - by default
                ### centers of clouds
            # color_g_mean = '#795fd7'
            color_g_mean = 'midnightblue'
            color_e_mean = 'maroon'
            color_dist = 'gold'
            vector_bw_blobs = 0.7
                ### vectors from void_point to centers
            color_g_vector = 'gold'
            color_e_vector ='gold'
            vector_state_lw = 0.7
                ### zero points
            color_void = 'gold'
            color_zero = 'k'
            color_zero_vector = 'k'
            vector_zero_lw = 0.5
                ### background of image
            fig_face_color = 'white' #this does not work
            fig_border_color = 'r'
            bg_color = 'white'
            grid_color =  my_colors_dict['meduza_dark']
            grid_transp = 0.5
            title_color = 'k'
            legend_color = 'white'
            legend_text_color = 'k'
            legend_alpha = 0.7
            legend_frame_color = 'k'

            th_alpha = 0.7
            th_color = my_colors_dict['deus_ex_gold']
            th_width = 2.0
            th_linestyle = '--'


            import matplotlib as mpl
            AXES_COLOR = my_colors_dict['meduza_dark']
            mpl.rc('axes', edgecolor=AXES_COLOR, labelcolor=AXES_COLOR, grid=True)
            mpl.rc('xtick', color=AXES_COLOR)
            mpl.rc('ytick', color=AXES_COLOR)
            mpl.rc('grid', color=AXES_COLOR)

        ### font default
        if font is None:
            plt.rc('font', family = 'Verdana')

        ###############################################
        ### all data plot here ______________________##
        ###############################################
        if figsize is None:
            fig, ax = plt.subplots(1, 1, sharey=True, tight_layout=True, facecolor=fig_face_color, edgecolor = fig_border_color)
        else:
            if (type(figsize) != list):
                print 'figsize must be a lsit'
                return False
            else:
                if len(figsize) != 2:
                    print 'figsize list must contain to numbers (x and y size)'
                    return False
            fig, ax = plt.subplots(1, 1, figsize=(figsize[0],figsize[1]), sharey=True, tight_layout=True, facecolor=fig_face_color, edgecolor = fig_border_color)


        maxval = 1 ##this variable we use for define plt.ylim()

        if self.threshold_r is not None:
            plt.axvline(x=self.threshold_r, alpha=th_alpha, c=th_color, lw=th_width, ls=th_linestyle)

        ### fake plot just for string in legend
        str_params = ''
        if self.dict_param is not None:
            str_params = 'Rdt:'+ my_stround( self.dict_param['rudat'],4 )+'dB; t_read:'+my_stround( 1e9*self.dict_param['t_read'],3 )+'ns'
            plt.plot([],[], label = str_params, visible=False)

        S_e = self.S_eff_e
        S_g = self.S_eff_g
        str_S_e = '; Se='+ str( round(S_e,2) )
        str_S_g = '; Sg='+ str( round(S_g,2) )

        S_e_selected = self.S_eff_e_selected
        S_g_selected = self.S_eff_g_selected
        str_S_e_selected = '; Se_sel='+ str( round(S_e_selected,2) )
        str_S_g_selected = '; Sg_sel='+ str( round(S_g_selected,2) )

        if regime == 'raw_data':
            print 'regime: raw_data'
            if (self.center_r_g is not None) and (self.center_r_e is not None):
                plt.axvline(x=self.center_r_g, alpha=th_alpha, c=color_g_mean, lw=th_width, ls=th_linestyle)
                plt.axvline(x=self.center_r_e, alpha=th_alpha, c=color_e_mean, lw=th_width, ls=th_linestyle)
            #### plot r_g, r_e hists ####
            if (self.hist_r_g.hist_xy is not None) and (self.hist_r_e.hist_xy is not None):
                hist_g   = self.hist_r_g.hist_xy
                hist_e   = self.hist_r_e.hist_xy
                maxval = np.max([ maxval, np.max(hist_g[0]), np.max(hist_e[0]) ])
                plt.plot(hist_g[1], hist_g[0], drawstyle='steps', lw=1*lw, color=color_g, label='Read g-state'+str_S_g)
                plt.plot(hist_e[1], hist_e[0], drawstyle='steps', lw=1*lw, color=color_e, label='Read e-state'+str_S_e)

            #### plot fit hists ####
            if (self.hist_r_g.gauss_fit is not None) and (self.hist_r_e.gauss_fit is not None):
                hist_g  = self.hist_r_g.gauss_fit
                hist_e  = self.hist_r_e.gauss_fit
                # maxval = np.max([ maxval, np.max(hist_g_fit[0]), np.max(hist_e_fit[0]) ])
                plt.plot(hist_g[1], hist_g[0], lw=1*lw, color=color_fit_g)
                plt.plot(hist_e[1], hist_e[0], lw=1*lw, color=color_fit_e)

        elif regime == 'raw_data_and_pre':
            print 'regime: raw_data_and_pre'
            if (self.center_r_g is not None) and (self.center_r_e is not None):
                plt.axvline(x=self.center_r_g, alpha=th_alpha, c=color_g_mean, lw=th_width, ls=th_linestyle)
                plt.axvline(x=self.center_r_e, alpha=th_alpha, c=color_e_mean, lw=th_width, ls=th_linestyle)
            #### plot PREPULSE hists ####
            if (self.hist_r_g_pre is not None) and (self.hist_r_e_pre is not None):
                hist_g  = self.hist_r_g_pre.hist_xy
                hist_e  = self.hist_r_e_pre.hist_xy
                # maxval = np.max([ maxval, np.max(hist_g_fit[0]), np.max(hist_e_fit[0]) ])
                plt.plot(hist_g[1], hist_g[0], drawstyle='steps', lw=1*lw, color=color_post, label='Prepulse g-state')
                plt.plot(hist_e[1], hist_e[0], drawstyle='steps', lw=1*lw, color=color_post, label='Prepulse e-state')

            #### plot r_g, r_e hists ####
            if (self.hist_r_g.hist_xy is not None) and (self.hist_r_e.hist_xy is not None):
                hist_g   = self.hist_r_g.hist_xy
                hist_e   = self.hist_r_e.hist_xy
                maxval = np.max([ maxval, np.max(hist_g[0]), np.max(hist_e[0]) ])
                plt.plot(hist_g[1], hist_g[0], drawstyle='steps', lw=1*lw, color=color_g, label='Read g-state'+str_S_g)
                plt.plot(hist_e[1], hist_e[0], drawstyle='steps', lw=1*lw, color=color_e, label='Read e-state'+str_S_e)

            #### plot fit hists ####
            if (self.hist_r_g.gauss_fit is not None) and (self.hist_r_e.gauss_fit is not None):
                hist_g  = self.hist_r_g.gauss_fit
                hist_e  = self.hist_r_e.gauss_fit
                # maxval = np.max([ maxval, np.max(hist_g_fit[0]), np.max(hist_e_fit[0]) ])
                plt.plot(hist_g[1], hist_g[0], lw=1*lw, color=color_fit_g)
                plt.plot(hist_e[1], hist_e[0], lw=1*lw, color=color_fit_e)


        elif regime=='selected':
            if (self.center_r_g_select is not None) and (self.center_r_e_select is not None):
                plt.axvline(x=self.center_r_g_select, alpha=th_alpha, c=color_g_mean, lw=th_width, ls=th_linestyle)
                plt.axvline(x=self.center_r_e_select, alpha=th_alpha, c=color_e_mean, lw=th_width, ls=th_linestyle)

            print 'regime: selected'
            #### plot r_g, r_e hists  SELECTED ####
            if (self.hist_r_g_select.hist_xy is not None) and (self.hist_r_e_select.hist_xy is not None):
                hist_g   = self.hist_r_g_select.hist_xy
                hist_e   = self.hist_r_e_select.hist_xy
                maxval = np.max([ maxval, np.max(hist_g[0]), np.max(hist_e[0]) ])
                plt.plot(hist_g[1], hist_g[0], drawstyle='steps', lw=1*lw, color=color_g, label='Selected g-state'+str_S_g_selected)
                plt.plot(hist_e[1], hist_e[0], drawstyle='steps', lw=1*lw, color=color_e, label='Selected e-state'+str_S_e_selected)

            #### plot fit hists ####
            if (self.hist_r_g_select.gauss_fit is not None) and (self.hist_r_e_select.gauss_fit is not None):
                hist_g  = self.hist_r_g_select.gauss_fit
                hist_e  = self.hist_r_e_select.gauss_fit
                # maxval = np.max([ maxval, np.max(hist_g_fit[0]), np.max(hist_e_fit[0]) ])
                plt.plot(hist_g[1], hist_g[0], lw=1*lw, color=color_fit_g)
                plt.plot(hist_e[1], hist_e[0], lw=1*lw, color=color_fit_e)


        elif regime=='raw_and_selected':
            if (self.center_r_g_select is not None) and (self.center_r_e_select is not None):
                plt.axvline(x=self.center_r_g_select, alpha=th_alpha, c=color_g_mean, lw=th_width, ls=th_linestyle)
                plt.axvline(x=self.center_r_e_select, alpha=th_alpha, c=color_e_mean, lw=th_width, ls=th_linestyle)
            print 'regime==raw_and_selected'
            #### plot r_g, r_e hists  SELECTED ####
            if (self.hist_r_g_select.hist_xy is not None) and (self.hist_r_e_select.hist_xy is not None):
                hist_g   = self.hist_r_g_select.hist_xy
                hist_e   = self.hist_r_e_select.hist_xy
                maxval = np.max([ maxval, np.max(hist_g[0]), np.max(hist_e[0]) ])
                plt.plot(hist_g[1], hist_g[0], drawstyle='steps', lw=1*lw, color='b', label='Selected g-state'+str_S_g_selected)
                plt.plot(hist_e[1], hist_e[0], drawstyle='steps', lw=1*lw, color='r', label='Selected e-state'+str_S_e_selected)

            #### plot r_g, r_e hists ####
            if (self.hist_r_g.hist_xy is not None) and (self.hist_r_e.hist_xy is not None):
                hist_g   = self.hist_r_g.hist_xy
                hist_e   = self.hist_r_e.hist_xy
                maxval = np.max([ maxval, np.max(hist_g[0]), np.max(hist_e[0]) ])
                plt.plot(hist_g[1], hist_g[0], drawstyle='steps', lw=1*lw, color=color_g, label='Read g-state'+str_S_g)
                plt.plot(hist_e[1], hist_e[0], drawstyle='steps', lw=1*lw, color=color_e, label='Read e-state'+str_S_e)

            #### plot fit hists ####
            ## selected
            if (self.hist_r_g_select.gauss_fit is not None) and (self.hist_r_e_select.gauss_fit is not None):
                hist_g  = self.hist_r_g_select.gauss_fit
                hist_e  = self.hist_r_e_select.gauss_fit
                # maxval = np.max([ maxval, np.max(hist_g_fit[0]), np.max(hist_e_fit[0]) ])
                plt.plot(hist_g[1], hist_g[0], lw=1*lw, color='b')
                plt.plot(hist_e[1], hist_e[0], lw=1*lw, color='r')
            ### raw
            if (self.hist_r_g.gauss_fit is not None) and (self.hist_r_e.gauss_fit is not None):
                hist_g  = self.hist_r_g.gauss_fit
                hist_e  = self.hist_r_e.gauss_fit
                # maxval = np.max([ maxval, np.max(hist_g_fit[0]), np.max(hist_e_fit[0]) ])
                plt.plot(hist_g[1], hist_g[0], lw=1*lw, color=color_fit_g)
                plt.plot(hist_e[1], hist_e[0], lw=1*lw, color=color_fit_e)

        else:
            print 'plot_hist: wrong value for regime!'

        ###############################################


        ###############################################
        ### all work with axes here _________________##
        ###############################################
        ax.set_facecolor(bg_color)
        plt.grid(color=grid_color, alpha= grid_transp)

        [leftlim, rightlim] = limits
        if leftlim is None:
            leftlim = np.min([ np.min(hist_g[1]), np.min(hist_e[1]) ])
        if rightlim is None:
            rightlim = np.max([ np.max(hist_g[1]), np.max(hist_e[1]) ])
        plt.xlim(leftlim,rightlim)

        if log:
            plt.yscale('log')
            plt.ylim(ymin=1, ymax=maxval*1.5)
        else:
            plt.yscale('linear')
            plt.ylim(ymin=1, ymax=maxval*1.1)

        if title_str is '':
            title_str = self.timestamp

        # if self.CONVERT_TOMV:
        #     lab_units = '[mV]'
        # else:
        #     lab_units = '[V]'
        lab_units = '[unitless]'

        if font is not None:
            plt.title(title_str, color=title_color,fontproperties = font)
            plt.xlabel(lab_units, fontproperties = font)
            plt.ylabel('Counts',fontproperties = font)
            leg = plt.legend(fancybox=True, framealpha=legend_alpha, loc='upper left', facecolor=legend_color, edgecolor=legend_frame_color,prop=font)
            for label in ax.get_xticklabels():  #set font to each xtick
                label.set_fontproperties(font)
            for label in ax.get_yticklabels():  #set font to each xtick
                label.set_fontproperties(font)
        else:
            plt.title(title_str, color=title_color)
            plt.xlabel(lab_units)
            plt.ylabel('Counts')
            leg = plt.legend(fancybox=True, framealpha=legend_alpha, loc='upper left', facecolor=legend_color, edgecolor=legend_frame_color)

        for text in leg.get_texts():        #set color to legend text
            plt.setp(text, color = legend_text_color)

        if save:
            if savepath == '':
                savepath='savings\\'
            import os
            if not os.path.exists(savepath):
                os.makedirs(savepath)
            full_fname = savepath +'\\'+ fname + '.png'
            if fig_transp:
                plt.savefig(full_fname, transparent = True)
            else:
                plt.savefig(full_fname,facecolor=fig_face_color, edgecolor=fig_border_color)

        if show:
            plt.show()
            return fig
        else:
            plt.close()
            return True

    def plot_f_vs_threshold(self, xmin=None, xmax=None, ymin=None, ymax=None):
        '''
        function to show fidelity versus threshold
        '''
        def get_fro_vs_threshold(x_g, x_e, threshold):
            '''
            Simplest function to calculate only f_ro for given threshold value
            Calculate fidelity for given value of THRESHOLD
            Takes 1D arrays of normalized g and e results. ( re_g & re_e )
            return only value of f_ro = 1. - 0.5*(p_ge + p_eg)
            '''
            dict_count = get_count_states(x_g,x_e,threshold)
            p_ge = dict_count['p_ge']
            p_eg = dict_count['p_eg']
            f_ro = 1. - 0.5*(p_ge + p_eg)
            return f_ro

        if (self.x_g is None) or (self.x_e is None):
            self.make_norm_data_from_raw()
        x_g = self.x_g
        x_e = self.x_e

        nop = 200.
        [leftlim, rightlim, rabbish1, rabbish2 ] = crop_fluctuations(self.x_g, [0] , self.x_e, [0], 0, 0 )
        thr_vector = np.linspace(leftlim, rightlim, nop)

        fid_vector = []
        for thr in thr_vector:
            fid_vector.append( get_fro_vs_threshold(x_g, x_e, thr) )

        pic = plt.figure()

        if self.CONVERT_TOMV:
            lab_units = '[mV]'
        else:
            lab_units = '[V]'

        plt.plot(thr_vector, fid_vector, '.')
        plt.xlabel('threshold '+lab_units)
        plt.ylabel('Fidelity')


        ##postselected
        if (self.x_g_select is not None ) and (self.x_e_select is not None):
            fid_post_vector = []
            for thr in thr_vector:
                fid_post_vector.append( get_fro_vs_threshold(self.x_g_select, self.x_e_select, thr) )
            plt.plot(thr_vector, fid_post_vector, '.')

        return pic

################################################################################

import numpy as np
import pandas as pd
from pymc import *

def bayesian_linear_regression(x_data, y_data, A_true=pd.DataFrame(np.zeros((46,46)))):
    N = len(y_data)
    F = np.shape(x_data)[1]
    prior_mean_ab01 = A_true.values[:,0]
    prior_mean_ab02 = A_true.values[:,1]
    prior_mean_ab03 = A_true.values[:,2]
    prior_mean_ab04 = A_true.values[:,3]
    prior_mean_ab05 = A_true.values[:,4]
    prior_mean_ab06 = A_true.values[:,5]
    prior_mean_ab07 = A_true.values[:,6]
    prior_mean_ab08 = A_true.values[:,7]
    prior_mean_ab09 = A_true.values[:,8]
    prior_mean_ab10 = A_true.values[:,9]
    prior_mean_ab11 = A_true.values[:,10]
    prior_mean_ab12 = A_true.values[:,11]
    prior_mean_ab13 = A_true.values[:,12]
    prior_mean_ab14 = A_true.values[:,13]
    prior_mean_ab15 = A_true.values[:,14]
    prior_mean_ab16 = A_true.values[:,15]
    prior_mean_ab17 = A_true.values[:,16]
    prior_mean_ab18 = A_true.values[:,17]
    prior_mean_ab19 = A_true.values[:,18]
    prior_mean_ab20 = A_true.values[:,19]
    prior_mean_ab21 = A_true.values[:,20]
    prior_mean_ab22 = A_true.values[:,21]
    prior_mean_ab23 = A_true.values[:,22]
    prior_mean_ab24 = A_true.values[:,23]
    prior_mean_ab25 = A_true.values[:,24]
    prior_mean_ab26 = A_true.values[:,25]
    prior_mean_ab27 = A_true.values[:,26]
    prior_mean_ab28 = A_true.values[:,27]
    prior_mean_ab29 = A_true.values[:,28]
    prior_mean_ab30 = A_true.values[:,29]
    prior_mean_ab31 = A_true.values[:,30]
    prior_mean_ab32 = A_true.values[:,31]
    prior_mean_ab33 = A_true.values[:,32]
    prior_mean_ab34 = A_true.values[:,33]
    prior_mean_ab35 = A_true.values[:,34]
    prior_mean_ab36 = A_true.values[:,35]
    prior_mean_ab37 = A_true.values[:,36]
    prior_mean_ab38 = A_true.values[:,37]
    prior_mean_ab39 = A_true.values[:,38]
    prior_mean_ab40 = A_true.values[:,39]
    prior_mean_ab41 = A_true.values[:,40]
    prior_mean_ab42 = A_true.values[:,41]
    prior_mean_ab43 = A_true.values[:,42]
    prior_mean_ab44 = A_true.values[:,43]
    prior_mean_ab45 = A_true.values[:,44]
    prior_mean_ab46 = A_true.values[:,45]
    scalar = 100
    offset_min = -10000
    offset_max = 10000
    
    beta0_ab01 = Uniform('offset_ab01', offset_min, offset_max)
    beta1_ab01 = MvNormal('weights_ab01', prior_mean_ab01, scalar*np.eye(F))
    beta0_ab02 = Uniform('offset_ab02', offset_min, offset_max)
    beta1_ab02 = MvNormal('weights_ab02', prior_mean_ab02, scalar*np.eye(F))
    beta0_ab03 = Uniform('offset_ab03', offset_min, offset_max)
    beta1_ab03 = MvNormal('weights_ab03', prior_mean_ab03, scalar*np.eye(F))
    beta0_ab04 = Uniform('offset_ab04', offset_min, offset_max)
    beta1_ab04 = MvNormal('weights_ab04', prior_mean_ab04, scalar*np.eye(F))
    beta0_ab05 = Uniform('offset_ab05', offset_min, offset_max)
    beta1_ab05 = MvNormal('weights_ab05', prior_mean_ab05, scalar*np.eye(F))
    beta0_ab06 = Uniform('offset_ab06', offset_min, offset_max)
    beta1_ab06 = MvNormal('weights_ab06', prior_mean_ab06, scalar*np.eye(F))
    beta0_ab07 = Uniform('offset_ab07', offset_min, offset_max)
    beta1_ab07 = MvNormal('weights_ab07', prior_mean_ab07, scalar*np.eye(F))
    beta0_ab08 = Uniform('offset_ab08', offset_min, offset_max)
    beta1_ab08 = MvNormal('weights_ab08', prior_mean_ab08, scalar*np.eye(F))
    beta0_ab09 = Uniform('offset_ab09', offset_min, offset_max)
    beta1_ab09 = MvNormal('weights_ab09', prior_mean_ab09, scalar*np.eye(F))
    beta0_ab10 = Uniform('offset_ab10', offset_min, offset_max)
    beta1_ab10 = MvNormal('weights_ab10', prior_mean_ab10, scalar*np.eye(F))
    beta0_ab11 = Uniform('offset_ab11', offset_min, offset_max)
    beta1_ab11 = MvNormal('weights_ab11', prior_mean_ab11, scalar*np.eye(F))
    beta0_ab12 = Uniform('offset_ab12', offset_min, offset_max)
    beta1_ab12 = MvNormal('weights_ab12', prior_mean_ab12, scalar*np.eye(F))
    beta0_ab13 = Uniform('offset_ab13', offset_min, offset_max)
    beta1_ab13 = MvNormal('weights_ab13', prior_mean_ab13, scalar*np.eye(F))
    beta0_ab14 = Uniform('offset_ab14', offset_min, offset_max)
    beta1_ab14 = MvNormal('weights_ab14', prior_mean_ab14, scalar*np.eye(F))
    beta0_ab15 = Uniform('offset_ab15', offset_min, offset_max)
    beta1_ab15 = MvNormal('weights_ab15', prior_mean_ab15, scalar*np.eye(F))
    beta0_ab16 = Uniform('offset_ab16', offset_min, offset_max)
    beta1_ab16 = MvNormal('weights_ab16', prior_mean_ab16, scalar*np.eye(F))
    beta0_ab17 = Uniform('offset_ab17', offset_min, offset_max)
    beta1_ab17 = MvNormal('weights_ab17', prior_mean_ab17, scalar*np.eye(F))
    beta0_ab18 = Uniform('offset_ab18', offset_min, offset_max)
    beta1_ab18 = MvNormal('weights_ab18', prior_mean_ab18, scalar*np.eye(F))
    beta0_ab19 = Uniform('offset_ab19', offset_min, offset_max)
    beta1_ab19 = MvNormal('weights_ab19', prior_mean_ab19, scalar*np.eye(F))
    beta0_ab20 = Uniform('offset_ab20', offset_min, offset_max)
    beta1_ab20 = MvNormal('weights_ab20', prior_mean_ab20, scalar*np.eye(F))
    beta0_ab21 = Uniform('offset_ab21', offset_min, offset_max)
    beta1_ab21 = MvNormal('weights_ab21', prior_mean_ab21, scalar*np.eye(F))
    beta0_ab22 = Uniform('offset_ab22', offset_min, offset_max)
    beta1_ab22 = MvNormal('weights_ab22', prior_mean_ab22, scalar*np.eye(F))
    beta0_ab23 = Uniform('offset_ab23', offset_min, offset_max)
    beta1_ab23 = MvNormal('weights_ab23', prior_mean_ab23, scalar*np.eye(F))
    beta0_ab24 = Uniform('offset_ab24', offset_min, offset_max)
    beta1_ab24 = MvNormal('weights_ab24', prior_mean_ab24, scalar*np.eye(F))
    beta0_ab25 = Uniform('offset_ab25', offset_min, offset_max)
    beta1_ab25 = MvNormal('weights_ab25', prior_mean_ab25, scalar*np.eye(F))
    beta0_ab26 = Uniform('offset_ab26', offset_min, offset_max)
    beta1_ab26 = MvNormal('weights_ab26', prior_mean_ab26, scalar*np.eye(F))
    beta0_ab27 = Uniform('offset_ab27', offset_min, offset_max)
    beta1_ab27 = MvNormal('weights_ab27', prior_mean_ab27, scalar*np.eye(F))
    beta0_ab28 = Uniform('offset_ab28', offset_min, offset_max)
    beta1_ab28 = MvNormal('weights_ab28', prior_mean_ab28, scalar*np.eye(F))
    beta0_ab29 = Uniform('offset_ab29', offset_min, offset_max)
    beta1_ab29 = MvNormal('weights_ab29', prior_mean_ab29, scalar*np.eye(F))
    beta0_ab30 = Uniform('offset_ab30', offset_min, offset_max)
    beta1_ab30 = MvNormal('weights_ab30', prior_mean_ab30, scalar*np.eye(F))
    beta0_ab31 = Uniform('offset_ab31', offset_min, offset_max)
    beta1_ab31 = MvNormal('weights_ab31', prior_mean_ab31, scalar*np.eye(F))
    beta0_ab32 = Uniform('offset_ab32', offset_min, offset_max)
    beta1_ab32 = MvNormal('weights_ab32', prior_mean_ab32, scalar*np.eye(F))
    beta0_ab33 = Uniform('offset_ab33', offset_min, offset_max)
    beta1_ab33 = MvNormal('weights_ab33', prior_mean_ab33, scalar*np.eye(F))
    beta0_ab34 = Uniform('offset_ab34', offset_min, offset_max)
    beta1_ab34 = MvNormal('weights_ab34', prior_mean_ab34, scalar*np.eye(F))
    beta0_ab35 = Uniform('offset_ab35', offset_min, offset_max)
    beta1_ab35 = MvNormal('weights_ab35', prior_mean_ab35, scalar*np.eye(F))
    beta0_ab36 = Uniform('offset_ab36', offset_min, offset_max)
    beta1_ab36 = MvNormal('weights_ab36', prior_mean_ab36, scalar*np.eye(F))
    beta0_ab37 = Uniform('offset_ab37', offset_min, offset_max)
    beta1_ab37 = MvNormal('weights_ab37', prior_mean_ab37, scalar*np.eye(F))
    beta0_ab38 = Uniform('offset_ab38', offset_min, offset_max)
    beta1_ab38 = MvNormal('weights_ab38', prior_mean_ab38, scalar*np.eye(F))
    beta0_ab39 = Uniform('offset_ab39', offset_min, offset_max)
    beta1_ab39 = MvNormal('weights_ab39', prior_mean_ab39, scalar*np.eye(F))
    beta0_ab40 = Uniform('offset_ab40', offset_min, offset_max)
    beta1_ab40 = MvNormal('weights_ab40', prior_mean_ab40, scalar*np.eye(F))
    beta0_ab41 = Uniform('offset_ab41', offset_min, offset_max)
    beta1_ab41 = MvNormal('weights_ab41', prior_mean_ab41, scalar*np.eye(F))
    beta0_ab42 = Uniform('offset_ab42', offset_min, offset_max)
    beta1_ab42 = MvNormal('weights_ab42', prior_mean_ab42, scalar*np.eye(F))
    beta0_ab43 = Uniform('offset_ab43', offset_min, offset_max)
    beta1_ab43 = MvNormal('weights_ab43', prior_mean_ab43, scalar*np.eye(F))
    beta0_ab44 = Uniform('offset_ab44', offset_min, offset_max)
    beta1_ab44 = MvNormal('weights_ab44', prior_mean_ab44, scalar*np.eye(F))
    beta0_ab45 = Uniform('offset_ab45', offset_min, offset_max)
    beta1_ab45 = MvNormal('weights_ab45', prior_mean_ab45, scalar*np.eye(F))
    beta0_ab46 = Uniform('offset_ab46', offset_min, offset_max)
    beta1_ab46 = MvNormal('weights_ab46', prior_mean_ab46, scalar*np.eye(F))

    @deterministic
    def mu_ab01(beta1_ab01=beta1_ab01, beta0_ab01=beta0_ab01, x_data=x_data):
        return np.dot(x_data, beta1_ab01.reshape((F,1))) + beta0_ab01
    @deterministic
    def mu_ab02(beta1_ab02=beta1_ab02, beta0_ab02=beta0_ab02, x_data=x_data):
        return np.dot(x_data, beta1_ab02.reshape((F,1))) + beta0_ab02
    @deterministic
    def mu_ab03(beta1_ab03=beta1_ab03, beta0_ab03=beta0_ab03, x_data=x_data):
        return np.dot(x_data, beta1_ab03.reshape((F,1))) + beta0_ab03
    @deterministic
    def mu_ab04(beta1_ab04=beta1_ab04, beta0_ab04=beta0_ab04, x_data=x_data):
        return np.dot(x_data, beta1_ab04.reshape((F,1))) + beta0_ab04
    @deterministic
    def mu_ab05(beta1_ab05=beta1_ab05, beta0_ab05=beta0_ab05, x_data=x_data):
        return np.dot(x_data, beta1_ab05.reshape((F,1))) + beta0_ab05
    @deterministic
    def mu_ab06(beta1_ab06=beta1_ab06, beta0_ab06=beta0_ab06, x_data=x_data):
        return np.dot(x_data, beta1_ab06.reshape((F,1))) + beta0_ab06
    @deterministic
    def mu_ab07(beta1_ab07=beta1_ab07, beta0_ab07=beta0_ab07, x_data=x_data):
        return np.dot(x_data, beta1_ab07.reshape((F,1))) + beta0_ab07
    @deterministic
    def mu_ab08(beta1_ab08=beta1_ab08, beta0_ab08=beta0_ab08, x_data=x_data):
        return np.dot(x_data, beta1_ab08.reshape((F,1))) + beta0_ab08
    @deterministic
    def mu_ab09(beta1_ab09=beta1_ab09, beta0_ab09=beta0_ab09, x_data=x_data):
        return np.dot(x_data, beta1_ab09.reshape((F,1))) + beta0_ab09
    @deterministic
    def mu_ab10(beta1_ab10=beta1_ab10, beta0_ab10=beta0_ab10, x_data=x_data):
        return np.dot(x_data, beta1_ab10.reshape((F,1))) + beta0_ab10
    @deterministic
    def mu_ab11(beta1_ab11=beta1_ab11, beta0_ab11=beta0_ab11, x_data=x_data):
        return np.dot(x_data, beta1_ab11.reshape((F,1))) + beta0_ab11
    @deterministic
    def mu_ab12(beta1_ab12=beta1_ab12, beta0_ab12=beta0_ab12, x_data=x_data):
        return np.dot(x_data, beta1_ab12.reshape((F,1))) + beta0_ab12
    @deterministic
    def mu_ab13(beta1_ab13=beta1_ab13, beta0_ab13=beta0_ab13, x_data=x_data):
        return np.dot(x_data, beta1_ab13.reshape((F,1))) + beta0_ab13
    @deterministic
    def mu_ab14(beta1_ab14=beta1_ab14, beta0_ab14=beta0_ab14, x_data=x_data):
        return np.dot(x_data, beta1_ab14.reshape((F,1))) + beta0_ab14
    @deterministic
    def mu_ab15(beta1_ab15=beta1_ab15, beta0_ab15=beta0_ab15, x_data=x_data):
        return np.dot(x_data, beta1_ab15.reshape((F,1))) + beta0_ab15
    @deterministic
    def mu_ab16(beta1_ab16=beta1_ab16, beta0_ab16=beta0_ab16, x_data=x_data):
        return np.dot(x_data, beta1_ab16.reshape((F,1))) + beta0_ab16
    @deterministic
    def mu_ab17(beta1_ab17=beta1_ab17, beta0_ab17=beta0_ab17, x_data=x_data):
        return np.dot(x_data, beta1_ab17.reshape((F,1))) + beta0_ab17
    @deterministic
    def mu_ab18(beta1_ab18=beta1_ab18, beta0_ab18=beta0_ab18, x_data=x_data):
        return np.dot(x_data, beta1_ab18.reshape((F,1))) + beta0_ab18
    @deterministic
    def mu_ab19(beta1_ab19=beta1_ab19, beta0_ab19=beta0_ab19, x_data=x_data):
        return np.dot(x_data, beta1_ab19.reshape((F,1))) + beta0_ab19
    @deterministic
    def mu_ab20(beta1_ab20=beta1_ab20, beta0_ab20=beta0_ab20, x_data=x_data):
        return np.dot(x_data, beta1_ab20.reshape((F,1))) + beta0_ab20
    @deterministic
    def mu_ab21(beta1_ab21=beta1_ab21, beta0_ab21=beta0_ab21, x_data=x_data):
        return np.dot(x_data, beta1_ab21.reshape((F,1))) + beta0_ab21
    @deterministic
    def mu_ab22(beta1_ab22=beta1_ab22, beta0_ab22=beta0_ab22, x_data=x_data):
        return np.dot(x_data, beta1_ab22.reshape((F,1))) + beta0_ab22
    @deterministic
    def mu_ab23(beta1_ab23=beta1_ab23, beta0_ab23=beta0_ab23, x_data=x_data):
        return np.dot(x_data, beta1_ab23.reshape((F,1))) + beta0_ab23
    @deterministic
    def mu_ab24(beta1_ab24=beta1_ab24, beta0_ab24=beta0_ab24, x_data=x_data):
        return np.dot(x_data, beta1_ab24.reshape((F,1))) + beta0_ab24
    @deterministic
    def mu_ab25(beta1_ab25=beta1_ab25, beta0_ab25=beta0_ab25, x_data=x_data):
        return np.dot(x_data, beta1_ab25.reshape((F,1))) + beta0_ab25
    @deterministic
    def mu_ab26(beta1_ab26=beta1_ab26, beta0_ab26=beta0_ab26, x_data=x_data):
        return np.dot(x_data, beta1_ab26.reshape((F,1))) + beta0_ab26
    @deterministic
    def mu_ab27(beta1_ab27=beta1_ab27, beta0_ab27=beta0_ab27, x_data=x_data):
        return np.dot(x_data, beta1_ab27.reshape((F,1))) + beta0_ab27
    @deterministic
    def mu_ab28(beta1_ab28=beta1_ab28, beta0_ab28=beta0_ab28, x_data=x_data):
        return np.dot(x_data, beta1_ab28.reshape((F,1))) + beta0_ab28
    @deterministic
    def mu_ab29(beta1_ab29=beta1_ab29, beta0_ab29=beta0_ab29, x_data=x_data):
        return np.dot(x_data, beta1_ab29.reshape((F,1))) + beta0_ab29
    @deterministic
    def mu_ab30(beta1_ab30=beta1_ab30, beta0_ab30=beta0_ab30, x_data=x_data):
        return np.dot(x_data, beta1_ab30.reshape((F,1))) + beta0_ab30
    @deterministic
    def mu_ab31(beta1_ab31=beta1_ab31, beta0_ab31=beta0_ab31, x_data=x_data):
        return np.dot(x_data, beta1_ab31.reshape((F,1))) + beta0_ab31
    @deterministic
    def mu_ab32(beta1_ab32=beta1_ab32, beta0_ab32=beta0_ab32, x_data=x_data):
        return np.dot(x_data, beta1_ab32.reshape((F,1))) + beta0_ab32
    @deterministic
    def mu_ab33(beta1_ab33=beta1_ab33, beta0_ab33=beta0_ab33, x_data=x_data):
        return np.dot(x_data, beta1_ab33.reshape((F,1))) + beta0_ab33
    @deterministic
    def mu_ab34(beta1_ab34=beta1_ab34, beta0_ab34=beta0_ab34, x_data=x_data):
        return np.dot(x_data, beta1_ab34.reshape((F,1))) + beta0_ab34
    @deterministic
    def mu_ab35(beta1_ab35=beta1_ab35, beta0_ab35=beta0_ab35, x_data=x_data):
        return np.dot(x_data, beta1_ab35.reshape((F,1))) + beta0_ab35
    @deterministic
    def mu_ab36(beta1_ab36=beta1_ab36, beta0_ab36=beta0_ab36, x_data=x_data):
        return np.dot(x_data, beta1_ab36.reshape((F,1))) + beta0_ab36
    @deterministic
    def mu_ab37(beta1_ab37=beta1_ab37, beta0_ab37=beta0_ab37, x_data=x_data):
        return np.dot(x_data, beta1_ab37.reshape((F,1))) + beta0_ab37
    @deterministic
    def mu_ab38(beta1_ab38=beta1_ab38, beta0_ab38=beta0_ab38, x_data=x_data):
        return np.dot(x_data, beta1_ab38.reshape((F,1))) + beta0_ab38
    @deterministic
    def mu_ab39(beta1_ab39=beta1_ab39, beta0_ab39=beta0_ab39, x_data=x_data):
        return np.dot(x_data, beta1_ab39.reshape((F,1))) + beta0_ab39
    @deterministic
    def mu_ab40(beta1_ab40=beta1_ab40, beta0_ab40=beta0_ab40, x_data=x_data):
        return np.dot(x_data, beta1_ab40.reshape((F,1))) + beta0_ab40
    @deterministic
    def mu_ab41(beta1_ab41=beta1_ab41, beta0_ab41=beta0_ab41, x_data=x_data):
        return np.dot(x_data, beta1_ab41.reshape((F,1))) + beta0_ab41
    @deterministic
    def mu_ab42(beta1_ab42=beta1_ab42, beta0_ab42=beta0_ab42, x_data=x_data):
        return np.dot(x_data, beta1_ab42.reshape((F,1))) + beta0_ab42
    @deterministic
    def mu_ab43(beta1_ab43=beta1_ab43, beta0_ab43=beta0_ab43, x_data=x_data):
        return np.dot(x_data, beta1_ab43.reshape((F,1))) + beta0_ab43
    @deterministic
    def mu_ab44(beta1_ab44=beta1_ab44, beta0_ab44=beta0_ab44, x_data=x_data):
        return np.dot(x_data, beta1_ab44.reshape((F,1))) + beta0_ab44
    @deterministic
    def mu_ab45(beta1_ab45=beta1_ab45, beta0_ab45=beta0_ab45, x_data=x_data):
        return np.dot(x_data, beta1_ab45.reshape((F,1))) + beta0_ab45
    @deterministic
    def mu_ab46(beta1_ab46=beta1_ab46, beta0_ab46=beta0_ab46, x_data=x_data):
        return np.dot(x_data, beta1_ab46.reshape((F,1))) + beta0_ab46

    ab01 = Normal('ab01', mu=mu_ab01, tau=1, value=y_data[:,0], observed=True)
    ab02 = Normal('ab02', mu=mu_ab02, tau=1, value=y_data[:,1], observed=True)
    ab03 = Normal('ab03', mu=mu_ab03, tau=1, value=y_data[:,2], observed=True)
    ab04 = Normal('ab04', mu=mu_ab04, tau=1, value=y_data[:,3], observed=True)
    ab05 = Normal('ab05', mu=mu_ab05, tau=1, value=y_data[:,4], observed=True)
    ab06 = Normal('ab06', mu=mu_ab06, tau=1, value=y_data[:,5], observed=True)
    ab07 = Normal('ab07', mu=mu_ab07, tau=1, value=y_data[:,6], observed=True)
    ab08 = Normal('ab08', mu=mu_ab08, tau=1, value=y_data[:,7], observed=True)
    ab09 = Normal('ab09', mu=mu_ab09, tau=1, value=y_data[:,8], observed=True)
    ab10 = Normal('ab10', mu=mu_ab10, tau=1, value=y_data[:,9], observed=True)
    ab11 = Normal('ab11', mu=mu_ab11, tau=1, value=y_data[:,10], observed=True)
    ab12 = Normal('ab12', mu=mu_ab12, tau=1, value=y_data[:,11], observed=True)
    ab13 = Normal('ab13', mu=mu_ab13, tau=1, value=y_data[:,12], observed=True)
    ab14 = Normal('ab14', mu=mu_ab14, tau=1, value=y_data[:,13], observed=True)
    ab15 = Normal('ab15', mu=mu_ab15, tau=1, value=y_data[:,14], observed=True)
    ab16 = Normal('ab16', mu=mu_ab16, tau=1, value=y_data[:,15], observed=True)
    ab17 = Normal('ab17', mu=mu_ab17, tau=1, value=y_data[:,16], observed=True)
    ab18 = Normal('ab18', mu=mu_ab18, tau=1, value=y_data[:,17], observed=True)
    ab19 = Normal('ab19', mu=mu_ab19, tau=1, value=y_data[:,18], observed=True)
    ab20 = Normal('ab20', mu=mu_ab20, tau=1, value=y_data[:,19], observed=True)
    ab21 = Normal('ab21', mu=mu_ab21, tau=1, value=y_data[:,19], observed=True)
    ab22 = Normal('ab22', mu=mu_ab22, tau=1, value=y_data[:,19], observed=True)
    ab23 = Normal('ab23', mu=mu_ab23, tau=1, value=y_data[:,19], observed=True)
    ab24 = Normal('ab24', mu=mu_ab24, tau=1, value=y_data[:,19], observed=True)
    ab25 = Normal('ab25', mu=mu_ab25, tau=1, value=y_data[:,19], observed=True)
    ab26 = Normal('ab26', mu=mu_ab26, tau=1, value=y_data[:,19], observed=True)
    ab27 = Normal('ab27', mu=mu_ab27, tau=1, value=y_data[:,19], observed=True)
    ab28 = Normal('ab28', mu=mu_ab28, tau=1, value=y_data[:,19], observed=True)
    ab29 = Normal('ab29', mu=mu_ab29, tau=1, value=y_data[:,19], observed=True)
    ab30 = Normal('ab30', mu=mu_ab30, tau=1, value=y_data[:,19], observed=True)
    ab31 = Normal('ab31', mu=mu_ab31, tau=1, value=y_data[:,19], observed=True)
    ab32 = Normal('ab32', mu=mu_ab32, tau=1, value=y_data[:,19], observed=True)
    ab33 = Normal('ab33', mu=mu_ab33, tau=1, value=y_data[:,19], observed=True)
    ab34 = Normal('ab34', mu=mu_ab34, tau=1, value=y_data[:,19], observed=True)
    ab35 = Normal('ab35', mu=mu_ab35, tau=1, value=y_data[:,19], observed=True)
    ab36 = Normal('ab36', mu=mu_ab36, tau=1, value=y_data[:,19], observed=True)
    ab37 = Normal('ab37', mu=mu_ab37, tau=1, value=y_data[:,19], observed=True)
    ab38 = Normal('ab38', mu=mu_ab38, tau=1, value=y_data[:,19], observed=True)
    ab39 = Normal('ab39', mu=mu_ab39, tau=1, value=y_data[:,19], observed=True)
    ab40 = Normal('ab40', mu=mu_ab40, tau=1, value=y_data[:,19], observed=True)
    ab41 = Normal('ab41', mu=mu_ab41, tau=1, value=y_data[:,19], observed=True)
    ab42 = Normal('ab42', mu=mu_ab42, tau=1, value=y_data[:,19], observed=True)
    ab43 = Normal('ab43', mu=mu_ab43, tau=1, value=y_data[:,19], observed=True)
    ab44 = Normal('ab44', mu=mu_ab44, tau=1, value=y_data[:,19], observed=True)
    ab45 = Normal('ab45', mu=mu_ab45, tau=1, value=y_data[:,19], observed=True)
    ab46 = Normal('ab46', mu=mu_ab46, tau=1, value=y_data[:,19], observed=True)
    
    return locals()

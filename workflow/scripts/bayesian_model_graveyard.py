# using this as a graveyard for old models that I have updated away

# 220325 retiring first one, as it couldn't deal with >1 nucleosome
# initialize model
with Model() as smf:
    # Hyperparameters
    p_unmeth_given_bound = 0.99 # 0.9
    p_unmeth_given_unbound = 0.01 # 0.1
    # t is the thing we measure as a 1 in the data
    p_t_given_meth = 0.1
    p_t_given_unmeth = 0.99
    
    # marginalize
    p_t_given_bound = p_unmeth_given_bound * p_t_given_unmeth + (1 - p_unmeth_given_bound) * p_t_given_meth
    p_t_given_unbound = p_unmeth_given_unbound * p_t_given_unmeth + (1 - p_unmeth_given_unbound) * p_t_given_meth
        
    # nucleosome
    nuc_present = Bernoulli('nuc_present', p=0.5)
    nuc_start = DiscreteUniform('nuc_start', lower=200, upper=490)
    nuc_length = Normal('nuc_length', mu=147, sigma=7)
    
    # tfbs
    p_bound = 0.6 # can also modify this to have site-specific prior, EB on the dataset?
    
    # do conditional on nucleosome position to get p_tf_bound_i
    binding_probs = [switch(and_(nuc_present>0.5,and_(nuc_start<=tfbs_positions[idx][1],nuc_start+nuc_length>=[idx][0])),0.01,p_bound) for idx in range(n_tfbs)]

    # make a bernoulli RV from p_bound
    binding_bernoullis = [Bernoulli('b{}'.format(idx+1), p=binding_probs[idx]) for idx in range(n_tfbs)]
    
    # make masks
    masks = [[1 if (p>=tfbs_positions[idx][0] and p<=tfbs_positions[idx][1]) else 0 for p in methyl_positions] for idx in range(n_tfbs)]
    series_masks = [pd.Series(masks[idx], index=methyl_positions) for idx in range(n_tfbs)]

    # convert the states into occupancies at each site
    tf_bound = 0
    for idx in range(n_tfbs):
        tf_bound += binding_bernoullis[idx] * series_masks[idx]

    nuc_bound = switch(and_(nuc_present>0.5,and_(methyl_positions>=nuc_start,methyl_positions<=nuc_start+nuc_length)),1,0)
    bound_bernoulli_p = switch(or_(tf_bound,nuc_bound),p_t_given_bound,p_t_given_unbound)

    # data
    data = Data("data", observed_data[0])
    
    # likelihood
    observed_seq = Bernoulli('observed_seq', p=bound_bernoulli_p, observed=data)

def convert_posterior_to_assignment(trace, ntfbs, nuc_thresh=0.5, tfbs_thresh=0.4):
    classification = {}

    nuc_present = np.array(trace.posterior.nuc_present).mean() > nuc_thresh
    classification['nuc_present'] = nuc_present

    nuc_start = int(np.array(trace.posterior.nuc_start).mean())
    nuc_end = nuc_start + int(np.array(trace.posterior.nuc_length).mean())
    classification['nuc_start'] = nuc_start
    classification['nuc_end'] = nuc_end

    for b in range(1,ntfbs+1):
        tfbs_state = np.array(getattr(trace.posterior, 'b{}'.format(b))).mean() > tfbs_thresh
        classification['tfbs_{}'.format(b)] = tfbs_state

    return classification
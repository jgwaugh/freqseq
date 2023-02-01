def get_p_sucess(p: float, delta: float) -> float:
    """
    Gets the probaiblity of stepping up given a relative
    difference between two rates

    Parameters
    ----------
    p: float
        Probability of being assigned to the treatment group
    delta: float
        The relative difference between treatment and control
        p_t = (1 + delta)p_c

    Returns
    -------
    float
        p(treatment = 1 | treatment + control = 1)

    """

    return p * (1 + delta) / (1 + p * delta)

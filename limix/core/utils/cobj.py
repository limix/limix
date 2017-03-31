class cObject(object):
    """
    This class add to object a util function for caching
    """
    def __init__(self):
        pass

    def clear_cache(self,*method_names):
        for method_name in method_names:
            setattr(self,'_cache_'+method_name,None)
            setattr(self,'_cached_'+method_name,False)

    def fill_cache(self,method_name,value):
        setattr(self,'_cache_'+method_name,value)
        setattr(self,'_cached_'+method_name,True)

    def clear_cache_idxs(self,method_name,idxs):
        if type(idxs) is not tuple: idxs = [idxs]
        tail = '_'.join(str(idx) for idx in idxs)
        setattr(self,'_cache_'+method_name+'_'+tail,None)
        setattr(self,'_cached_'+method_name+'_'+tail,False)

    def fill_cache_idxs(self,method_name,idxs,value):
        if type(idxs) is not list: idxs = [idxs]
        tail = '_'.join(str(idx) for idx in idxs)
        setattr(self,'_cache_'+method_name+'_'+tail,value)
        setattr(self,'_cached_'+method_name+'_'+tail,True)

def cached(method):
    """ this function is used as a decorator for caching """
    _cache_attr_name = '_cache_'+method.__name__
    _bool_attr_name  = '_cached_'+method.__name__
    def method_wrapper(self,*args,**kwargs):
        is_cached = getattr(self,_bool_attr_name)
        if not is_cached:
            result = method(self, *args, **kwargs)
            setattr(self, _cache_attr_name, result)
            setattr(self, _bool_attr_name, True)
        return getattr(self,'_cache_'+method.__name__)
    return method_wrapper

def cached_idxs(method):
    """ this function is used as a decorator for caching """
    def method_wrapper(self,*args,**kwargs):
        tail = '_'.join(str(idx) for idx in args)
        _cache_attr_name = '_cache_'+method.__name__+'_'+tail
        _bool_attr_name  = '_cached_'+method.__name__+'_'+tail
        is_cached = getattr(self,_bool_attr_name)
        if not is_cached:
            result = method(self, *args, **kwargs)
            setattr(self, _cache_attr_name, result)
            setattr(self, _bool_attr_name, True)
        return getattr(self,_cache_attr_name)
    return method_wrapper


'''
Author: Rui Qin
Date: 2025-09-29 17:06:01
LastEditTime: 2025-09-29 17:18:27
Description: 
'''
import operator

ops = {
    '>=': operator.ge,
    '<=': operator.le,
    '>': operator.gt,
    '<': operator.lt,
    '==': operator.eq,
    '!=': operator.ne
}

def filter(value, condition_dict) -> bool:
    if value is None:
        return False
    operator, thres = condition_dict
    return ops[operator](value, thres)

def screen_with_stats(results:list, conditions:dict) -> tuple[list, dict]:
    """Screen evalution results with conditions and return filter stats.
    """
    filtered = results.copy()
    filter_stats = {}
    
    for cls, subclses in conditions.items():
        filter_stats[cls] = {}
        for metrics in subclses.values():
            for m in metrics:
                filter_stats[cls][m] = 0
    
    for cls, subclses in conditions.items():
        for subcls, metrics in subclses.items():
            for metric, condition in metrics.items():
                new_filtered = []
                
                for result in filtered:
                    if cls not in result or result[cls] is None:
                        new_filtered = filtered
                        filter_stats[cls][metric] = len(new_filtered)
                        continue
                    if subcls not in result[cls]:
                        continue
                        
                    cls_data = result[cls][subcls]
                    if metric not in cls_data:
                        raise ValueError(f"Metric '{metric}' not found in results \
                                         for class '{cls}' and subclass '{subcls}'.")
                    
                    value = cls_data[metric]
                    
                    if filter(value, condition):
                        new_filtered.append(result)
                
                filter_stats[cls][metric] = len(new_filtered)
                filtered = new_filtered
    return filtered, filter_stats
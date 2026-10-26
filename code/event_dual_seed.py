"""Translate earlier point-identity dual penalties into threshold weights.

This only provides search weights. It does not promote the old certificate or
any new pricing bounds to an independently checked proof.
"""


def convert(source,certificate,scale=10**6):
    R,n = source['core_residual'],len(source['high_primes'])
    limits = [14 if h<n else 16 for h in range((1<<n)-1)]
    arrays = [[[0]*(limit+1) for x in R] for limit in limits]
    denominator = certificate['denominator']
    def coefficient(y):
        assert type(y) is int and y<=0
        return y*scale//denominator
    for x,line in enumerate(certificate['link_multipliers']):
        for h,y in enumerate(line):
            for t in range(1,limits[h]+1):
                arrays[h][x][t] += coefficient(y)
    for row in certificate['tail_multipliers']:
        h,y = row['axis'],coefficient(row['multiplier'])
        for x,v in enumerate(R):
            for t in range(row['thresholds'][v%2]+1,limits[h]+1):
                arrays[h][x][t] += y
    for row in certificate.get('rank_multipliers',[]):
        h,y = row['axis'],coefficient(row['multiplier'])
        for x,v in enumerate(R):
            arrays[h][x][row['thresholds'][v%2]] += y*row['weights'][v%2]
    for kind in ('pair_multipliers','subset_multipliers'):
        for row in certificate.get(kind,[]):
            h,y = row['axis'],coefficient(row['multiplier'])
            for x,t in zip(row['points'],row['thresholds']):
                arrays[h][x][t] += y
    for row in certificate.get('group_event_multipliers',[]):
        h,y = row['axis'],coefficient(row['multiplier'])
        for event in row['events']:
            arrays[h][event['point']][event['threshold']] += y*event['weight']
    return arrays

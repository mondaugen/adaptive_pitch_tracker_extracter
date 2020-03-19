from attack_finder import event_closeness_limiter

events=[1,6,2,3,9,1,5,4,4]
print(sorted(events))
print(event_closeness_limiter(events,2))

events2=[1,1,1,1,1,1,1]
print(event_closeness_limiter(events2,2))

events3=[1,2,3,4,5,6,7]
print(event_closeness_limiter(events3,2))

events4=[1]
print(event_closeness_limiter(events4,2))

events5=[]
print(event_closeness_limiter(events5,2))

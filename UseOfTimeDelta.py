# Basic timedelta usage for finance
import datetime as dt
from datetime import timedelta

# generic dates as input
initialdate = dt.date(2024, 6, 18)
enddate = dt.date(2025, 7, 18)

# difference between them as a timedelta
deltatime = enddate - initialdate

# adjust according to the required basis
deltatime_Act365 = deltatime.days / 365.25
deltatime_Act360 = deltatime.days / 360



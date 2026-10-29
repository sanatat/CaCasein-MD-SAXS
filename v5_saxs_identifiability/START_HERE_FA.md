# راهنمای فارسی CaCasein-MD-SAXS نسخه ۵

## هدف

هدف نسخه ۵ این است که نشان دهیم چرا SAXS نمی‌تواند همیشه ساختار واقعی را به‌صورت یکتا تعیین کند و چطور می‌توان با مدل‌سازی و گزارش عدم‌قطعیت نتیجه بهتری ارائه داد.

**نکته صداقت علمی:** در این نسخه شبیه‌سازی MD جدید اجرا نشده است. ۱۵ شبیه‌سازی واقعی نسخه ۴ دوباره تحلیل شده‌اند و آزمون‌های جدید معادله Debye، حل مسئله معکوس و کره تحلیلی انجام شده‌اند. Shape2SAS اصلی و داده‌های SAXS تجربی در این اجرا در دسترس نبودند؛ ما ادعای انجام آنها را نداریم.

## اجرا روی لینوکس

```bash
unzip CaCasein-MD-SAXS-v5-complete.zip
cd CaCasein-MD-SAXS/v5_saxs_identifiability
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements_v5.txt
python -m pytest -q tests
python run_v5.py --v4 ../v4_validation --out results_recomputed
```

برای اجرای همه ۳۶ تست پروژه که در محیط حاضر موفق بودند:

```bash
cd ..
python -m pytest -q tests v2_sequence_calcium/tests v3_explicit_ca/tests v4_validation/tests v5_saxs_identifiability/tests
```

## نتایج اصلی

- اعتبارسنجی Debye در برابر فرمول تحلیلی کره: خطای RMS مطلق ۰٫۰۰۲۹۴۵ برای شدت نرمال‌شده.
- دو ساختار متفاوت A و B، پانزده فاصله دوتایی یکسان و منحنی SAXS برابر با اختلاف حداکثر ۳٫۳۳×۱۰⁻۱۶ دارند.
- در ۶۰ تکرار داده‌سازی با نویز و احتمال‌های اولیه برابر، احتمال مدل A و B هرکدام تقریباً ۵۰٪ است.
- نتایج وابستگی به شرایط اولیه اتصال کلسیم از نسخه ۴ همچنان مشاهده می‌شود؛ مدل کلسیم تأیید تجربی ندارد.

## چه فایل‌هایی مهم هستند؟

- `run_v5.py`: محاسبات اصلی
- `saxs_core.py`: معادله Debye، P(r)، فرمول کره، Rg، Guinier و احتمال شرطی مدل‌ها
- `results/01_sphere_benchmark.png`: اعتبارسنجی عددی
- `results/02_exact_homometry.png`: اثبات ابهام ساختار
- `results/03_inverse_uncertainty.png`: حل مسئله معکوس با نویز
- `results/04_v4_reanalysis.png`: بازتحلیل MD قبلی
- `results/results_summary.json`: اعداد دقیق نتایج
- `README.md`: گزارش انگلیسی، محدودیت‌ها و منابع
- `analyze_experimental_saxs.py`: تحلیل مقدماتی فایل واقعی *در صورت تهیه داده معتبر*
- `check_shape2sas_reference.py`: مقایسه اختیاری با داده خروجی نرم‌افزار رسمی Shape2SAS (خود آن را اجرا نکرده‌ایم)

## مناسب رزومه

می‌توانی بنویسی:

> Built and tested a reproducible SAXS modelling and inverse-inference pipeline using numerical Debye scattering, analytical sphere benchmarks, synthetic-noise ensembles, exact homometric structure examples, and retrospective analysis of reduced-unit coarse-grained MD trajectories; explicitly evaluated structural non-identifiability and sampling sensitivity.

ننویس که Martini 3 یا تحلیل SAXS تجربی یا اتصال واقعی کلسیم به پروتئین اعتبارسنجی شده؛ هنوز این کارها انجام نشده‌اند.

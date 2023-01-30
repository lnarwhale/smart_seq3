#!/bin/bash
source bulk.sh
diff="smart_mou_X2C_3-smart_mou_X0C|smart_za_X2C-smart_za_X0C|nss_za_X2C-nss_za_X0C"
tss_TPM_compare $diff /data1/shenluemou/ing_data/test/re /data1/shenluemou/ing_data/test/data /data1/shenluemou/ing_data/test/smart_seq/pwm mouse

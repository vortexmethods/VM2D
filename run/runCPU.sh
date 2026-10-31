#!/bin/bash
# Запуск CPU-версии VM2D на 14 потоках: P- и E-ядра (cpu0-13),
# без двух LP-ядер (cpu14-15), которые тормозят параллельные циклы.
# Использование: ./runCPU.sh [исполняемый файл, по умолчанию ./VM2D]
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-14}
exec taskset -c ${VM_CPUS:-0-13} ${1:-./VM2D}

Запуск VM2D и обработка результатов
===================================

Запуск VM2D
-----------
		  
Запуск программы осуществляется из *рабочего каталога*, т.е. каталога, содержащего файлы ``problems``, ``defaults``, ``switchers`` и ``mechanics``, а также папки конкретных решаемых задач. Для запуска расчетов необходимо вызвать команду 

      ``VM2D``

(при необходимости - с указанием полного пути к скомпилированному исполняемому файлу; этот путь можно указать в системной переменной ``PATH``)


Сохранение результатов
----------------------

Вихревые следы сохраняются в vtk-файлы в подкаталог ``snapshots``, поля скоростей и давления, вычисленные в заданных точках --- в подкаталог ``velPress``. Обе этих папки располагаются в подкаталоге задачи и создаются автоматически. 

Также в подкаталог задачи сохраняется информация о нагрузках, действующих на обтекаемый профиль или на каждый профиль из системы обтекаемых профилей. Эта информация сохраняется в файлы с именами ``forces-airfoil-n``, где ``n`` соответствует номеру профиля по порядку из их списка, перечисленного в ключе *airfoil* в *паспорте задачи*. Данные файлы сохраняются в двух форматах: текстовом (без расширения) и *csv*.

Если профиль подвижен, то информация о его положении сохраняется в файл ``position-airfoil-n`` с той же логикой именования файлов, также в двух форматах.

Для просмотра vtk и csv-файлов удобно использовать свободно распространяемый постпроцессор с открытым исходным кодом *paraview*.


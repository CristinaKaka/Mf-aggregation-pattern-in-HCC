/* 
* Script created by Yulan Weng
* rename the annotations using the name of the TMA core
*/
for (core in getTMACoreList()) {
   core.getChildObjects().each {it.setName(core.getName())}
}
fireHierarchyUpdate()